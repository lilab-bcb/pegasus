import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

import numba
from numba import njit

from typing import List, Optional, Union
from pegasusio import UnimodalData, MultimodalData
from pegasus.tools import eff_n_jobs

import anndata

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer

def _gen_anndata(
    data: Union[MultimodalData, UnimodalData],
    features: str = "highly_variable_features",
    obs_columns: Optional[List[str]] = None,
    matkey: str = "raw.X",
) -> anndata.AnnData:
    """ Generate a new Anndata object for scvitools

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    features: ``str``, optional, default: ``"highly_variable_features"``
        Keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.
    obs_columns: ``List[str]``
        A list of obs keys that should be included in the new anndata.
    matkey: ``str``, optional, default: ``"raw.X"``
        Matrix key for the raw count

    Returns
    -------
    An AnnData object.
    """
    mat = data.get_matrix(matkey)
    if obs_columns is not None and obs_columns:
        obs_field = data.obs[obs_columns]
    else:
        obs_field = data.obs
    if features is not None:
        assert features in data.var
        idx = data.var[features].values
        X = mat[:, idx]
        var_field = data.var.loc[idx, []]
    else:
        X = mat
        var_field = data.var[[]]

    return anndata.AnnData(X = X, obs = obs_field, var = var_field)


@njit(fastmath=True, cache=True)
def _select_csr(data, indices, indptr, indexer, new_size):
    data_new = np.zeros_like(data[0:new_size])
    indices_new = np.zeros_like(indices[0:new_size])
    indptr_new = np.zeros_like(indptr)

    cnt = 0
    for i in range(indptr.size - 1):
        indptr_new[i] = cnt
        for j in range(indptr[i], indptr[i + 1]):
            new_idx = indexer[indices[j]]
            if new_idx >= 0:
                data_new[cnt] = data[j]
                indices_new[cnt] = new_idx
                cnt += 1
    indptr_new[indptr.size - 1] = cnt

    return data_new, indices_new, indptr_new


def _gen_query_anndata(
    data: Union[MultimodalData, UnimodalData],
    ref_features: pd.Index,
    obs_columns: Optional[List[str]] = None,
    matkey: str = "raw.X",
) -> anndata.AnnData:
    """ Generate a new query Anndata object for scvitools

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    ref_features: ``pd.Index``
        A pandas index of reference feature names
    obs_columns: ``List[str]``
        A list of obs keys that should be included in the new anndata.
    matkey: ``str``, optional, default: ``"raw.X"``
        Matrix key for the raw count

    Returns
    -------
    An AnnData object.
    """
    mat = data.get_matrix(matkey)
    if obs_columns is not None and obs_columns:
        obs_field = data.obs[obs_columns]
    else:
        obs_field = data.obs
    var_field = pd.DataFrame(index = ref_features)

    indexer = ref_features.get_indexer(data.var_names)
    new_size = (indexer[mat.indices]>=0).sum()
    data_new, indices_new, indptr_new = _select_csr(mat.data, mat.indices, mat.indptr, indexer, new_size)
    X = csr_matrix((data_new, indices_new, indptr_new), shape = (mat.shape[0], ref_features.size))
    X.sort_indices()

    return anndata.AnnData(X = X, obs = obs_field, var = var_field)


@timer(logger=logger)
def run_scvi(
    data: Union[MultimodalData, UnimodalData],
    features: str = "highly_variable_features",
    matkey: str = "raw.X",
    n_jobs: int = -1,
    random_state: int = 0,
    max_epochs: Union[int, None] = None,
    batch: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    use_gpu: Union[str, int, bool, None] = None,
) -> str:
    """Run scVI embedding.

    This is a wrapper of `scvitools <https://github.com/scverse/scvi-tools>`_ package.

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.
    features: ``str``, optional, default: ``"highly_variable_features"``
        Keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.
    matkey: ``str``, optional, default: ``"raw.X"``
        Matrix key for the raw count
    n_jobs : ``int``, optional, default: ``-1``.
        Number of threads to use. ``-1`` refers to using all physical CPU cores.
    random_state: ``int``, optional, default: ``0``.
        Seed for random number generator
    max_epochs: ``int | None``, optional, default: ``None``.
        Maximum number of training epochs. Defaults to np.min([round((20000 / n_cells) * 400), 400])
    batch: ``str``, optional, default: ``None``.
        If only one categorical covariate, the obs key representing batches that should be corrected for, default is ``None``.
    categorical_covariate_keys: ``List[str]``
        If multiple categorical covariates, a list of obs keys listing categorical covariates that should be corrected for, default is ``None``.
    continuous_covariate_keys: ``List[str]``
        A list of obs keys listing continuous covariates that should be corrected for, default is ``None``.
    use_gpu: ``str | int | bool | None``
        Use default GPU if available (if None or True), or index of GPU to use (if int), or name of GPU (if str, e.g., ``cuda:0``), or use CPU (if False).

    Returns
    -------
    out_rep: ``str``
        The keyword in ``data.obsm`` referring to the embedding calculated by integrative NMF algorithm. out_rep is always equal to "scVI"

    Update ``data.obsm``:
        * ``data.obsm['X_scVI']``: The embedding calculated by scVI.

    Examples
    --------
    >>> pg.run_scvi(data, batch="Channel")
    >>> pg.run_scvi(data, categorical_covariate_keys=["cell_source", "donor"], continuous_covariate_keys=["percent_mito", "percent_ribo"])
    """
    try:
        import scvi
    except ImportError as e:
        import sys
        logger.error(f"{e}\nscvi-tools needed! Try 'pip install scvi-tools'.")
        sys.exit(-1)

    logger.info("Start embedding using scVI.")

    obs_columns = []
    if batch is not None and batch:
        obs_columns.append(batch)
    if categorical_covariate_keys is not None and categorical_covariate_keys:
        obs_columns.extend(categorical_covariate_keys)
    if continuous_covariate_keys is not None and continuous_covariate_keys:
        obs_columns.extend(continuous_covariate_keys)

    adata = _gen_anndata(data, features, obs_columns, matkey) # gen AnnData

    scvi.settings.num_threads = eff_n_jobs(n_jobs) # set n_jobs
    scvi.settings.seed = random_state # set random_state, see [here](https://docs.scvi-tools.org/en/stable/_modules/scvi/_settings.html) for more details.

    if max_epochs is None:
        max_epochs = np.min([round((20000 / len(adata.obs)) * 400), 400])

    scvi.model.SCVI.setup_anndata(adata,
        batch_key=batch,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
    ) # register anndata

    model = scvi.model.SCVI(adata) # obtain scvi model
    model.train(max_epochs=max_epochs, use_gpu=use_gpu) # train model
    data.obsm["X_scVI"] = model.get_latent_representation() # Get embedding
    data.register_attr("X_scVI", "embedding") # Register X_scVI as an embedding

    return "scVI"


@timer(logger=logger)
def train_scarches_scanvi(
    data: Union[MultimodalData, UnimodalData],
    dir_path: str,
    label: str,
    unlabeled_category: str = "Unknown",
    features: str = "highly_variable_features",
    matkey: str = "raw.X",
    n_jobs: int = -1,
    random_state: int = 0,
    max_epochs: Union[int, None] = None,
    batch: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    semisupervised_max_epochs: Union[int, None] = None,
    n_samples_per_label: Optional[int] = None,
    use_gpu: Union[str, int, bool, None] = None,
    arches_params: dict = dict(
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
    ),
) -> None:
    """Run scArches training.

    This is a wrapper of `scvitools <https://github.com/scverse/scvi-tools>`_ package.

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.
    dir_path: ``str``.
        Save the model to this directory.
    label: ``str``.
        The obs key representing labels.
    unlabeled_category: ``str``, default: ``"Unknown"``
        Value used for unlabeled cells in ``label``.
    features: ``str``, optional, default: ``"highly_variable_features"``
        Keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.
    matkey: ``str``, optional, default: ``"raw.X"``
        Matrix key for the raw count
    n_jobs : ``int``, optional, default: ``-1``.
        Number of threads to use. ``-1`` refers to using all physical CPU cores.
    random_state: ``int``, optional, default: ``0``.
        Seed for random number generator
    max_epochs: ``int | None``, optional, default: ``None``.
        Maximum number of unsupervised training epochs. Defaults to np.min([round((20000 / n_cells) * 400), 400])
    batch: ``str``, optional, default: ``None``.
        If only one categorical covariate, the obs key representing batches that should be corrected for, default is ``None``.
    categorical_covariate_keys: ``List[str]``
        If multiple categorical covariates, a list of obs keys listing categorical covariates that should be corrected for, default is ``None``.
    continuous_covariate_keys: ``List[str]``
        A list of obs keys listing continuous covariates that should be corrected for, default is ``None``.
    semisupervised_max_epochs: ``int | None``, optional, default: ``None``.
        Maximum number of semisupervised training epochs. Defaults to np.min([round(np.sqrt(``max_epochs``)), 20])
    n_samples_per_label : ``int``, optional, default: ``None``.
        Number of subsamples for each label class to sample per epoch. By default, there is no label subsampling.
    use_gpu: ``str | int | bool | None``
        Use default GPU if available (if None or True), or index of GPU to use (if int), or name of GPU (if str, e.g., ‘cuda:0’), or use CPU (if False).
    arches_params: ``dict``.
        Hyperparameters for VAE. See https://docs.scvi-tools.org/en/stable/api/reference/scvi.module.VAE.html#scvi.module.VAE for more details

    Returns
    -------
    Update ``data.obsm``:
        * ``data.obsm['X_scVI']``: The embedding calculated by scVI.
        * ``data.obsm['X_scanVI']``: The embedding calculated by scanVI.

    Examples
    --------
    >>> pg.train_scarches_scanvi(data, dir_path="scanvi_model/", label="celltype", matkey="counts", batch="tech", n_samples_per_label=100)
    """
    try:
        import scvi
    except ImportError as e:
        import sys
        logger.error(f"{e}\nscvi-tools needed! Try 'pip install scvi-tools'.")
        sys.exit(-1)

    logger.info("Start training with scArches method.")

    obs_columns = [label]
    if batch is not None and batch:
        obs_columns.append(batch)
    if categorical_covariate_keys is not None and categorical_covariate_keys:
        obs_columns.extend(categorical_covariate_keys)
    if continuous_covariate_keys is not None and continuous_covariate_keys:
        obs_columns.extend(continuous_covariate_keys)

    adata = _gen_anndata(data, features, obs_columns, matkey) # gen AnnData

    scvi.settings.num_threads = eff_n_jobs(n_jobs) # set n_jobs
    scvi.settings.seed = random_state # set random_state, see [here](https://docs.scvi-tools.org/en/stable/_modules/scvi/_settings.html) for more details.

    # unsupervised
    if max_epochs is None:
        max_epochs = np.min([round((20000 / len(adata.obs)) * 400), 400])

    scvi.model.SCVI.setup_anndata(adata,
        batch_key=batch,
        categorical_covariate_keys=categorical_covariate_keys,
        continuous_covariate_keys=continuous_covariate_keys,
    ) # register anndata
    vae_ref = scvi.model.SCVI(adata, **arches_params) # obtain scvi model
    vae_ref.train(max_epochs=max_epochs, use_gpu=use_gpu) # train model

    data.obsm["X_scVI"] = vae_ref.get_latent_representation() # Get embedding
    data.register_attr("X_scVI", "embedding") # Register X_scVI as an embedding

    # semisupervised
    if semisupervised_max_epochs is None:
        semisupervised_max_epochs = np.min([round(np.sqrt(max_epochs)), 20])

    vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
        vae_ref,
        unlabeled_category=unlabeled_category,
        labels_key=label,
    )
    vae_ref_scan.train(
        max_epochs=semisupervised_max_epochs,
        n_samples_per_label=n_samples_per_label,
        use_gpu=use_gpu,
    )
    vae_ref_scan.save(dir_path, overwrite=True)

    data.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()
    data.register_attr("X_scANVI", "embedding")


@timer(logger=logger)
def predict_scarches_scanvi(
    data: Union[MultimodalData, UnimodalData],
    dir_path: str,
    label: str,
    predictions: str = "predictions",
    matkey: str = "raw.X",
    n_jobs: int = -1,
    random_state: int = 0,
    max_epochs: Union[int, None] = None,
    batch: Optional[str] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    use_gpu: Union[str, int, bool, None] = None,
) -> None:
    """Run scArches training.

    This is a wrapper of `scvitools <https://github.com/scverse/scvi-tools>`_ package.

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.
    dir_path: ``str``.
        Save the model to this directory.
    label: ``str``.
        The obs key representing labels.
    predictions: ``str``, , optional, default: ``"predictions"``
        The obs key to store predicted labels.
    matkey: ``str``, optional, default: ``"raw.X"``
        Matrix key for the raw count
    n_jobs : ``int``, optional, default: ``-1``.
        Number of threads to use. ``-1`` refers to using all physical CPU cores.
    random_state: ``int``, optional, default: ``0``.
        Seed for random number generator
    max_epochs: ``int | None``, optional, default: ``None``.
        Maximum number of training epochs. Defaults to np.min([round((20000 / n_cells) * 100), 100])
    batch: ``str``, optional, default: ``None``.
        If only one categorical covariate, the obs key representing batches that should be corrected for, default is ``None``.
    categorical_covariate_keys: ``List[str]``
        If multiple categorical covariates, a list of obs keys listing categorical covariates that should be corrected for, default is ``None``.
    continuous_covariate_keys: ``List[str]``
        A list of obs keys listing continuous covariates that should be corrected for, default is ``None``.
    use_gpu: ``str | int | bool | None``
        Use default GPU if available (if None or True), or index of GPU to use (if int), or name of GPU (if str, e.g., ‘cuda:0’), or use CPU (if False).

    Returns
    -------
    Update ``data.obsm``:
        * ``data.obsm['X_scanVI']``: The embedding calculated by scanVI.
        * ``data.obsm[predictions]``: The predicted labels by scanVI.

    Examples
    --------
    >>> pg.predict_scarches_scanvi(data, dir_path="scanvi_model/", label="celltype", matkey="counts", batch="tech")
    """
    try:
        import scvi
    except ImportError as e:
        import sys
        logger.error(f"{e}\nscvi-tools needed! Try 'pip install scvi-tools'.")
        sys.exit(-1)

    logger.info("Start prediction with scArches method.")

    obs_columns = [label]
    if batch is not None and batch:
        obs_columns.append(batch)
    if categorical_covariate_keys is not None and categorical_covariate_keys:
        obs_columns.extend(categorical_covariate_keys)
    if continuous_covariate_keys is not None and continuous_covariate_keys:
        obs_columns.extend(continuous_covariate_keys)

    features = scvi.model.SCANVI.prepare_query_anndata(None, dir_path, return_reference_var_names=True)
    adata = _gen_query_anndata(data, features, obs_columns, matkey) # gen AnnData

    scvi.settings.num_threads = eff_n_jobs(n_jobs) # set n_jobs
    scvi.settings.seed = random_state # set random_state, see [here](https://docs.scvi-tools.org/en/stable/_modules/scvi/_settings.html) for more details.

    if max_epochs is None:
        max_epochs = np.min([round((20000 / len(adata.obs)) * 100), 100])

    vae_q = scvi.model.SCANVI.load_query_data(adata, dir_path)
    vae_q.train(
        max_epochs=max_epochs,
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )

    data.obsm["X_scANVI"] = vae_q.get_latent_representation()
    data.register_attr("X_scANVI", "embedding")
    data.obs[predictions] = vae_q.predict()
