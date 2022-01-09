import time
import numpy as np
import pandas as pd
from typing import Union
from pegasusio import UnimodalData, MultimodalData

from pegasus.tools import select_features, X_from_rep, check_batch_key

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



@timer(logger=logger)
def run_harmony(
    data: Union[MultimodalData, UnimodalData],
    batch: str = "Channel",
    rep: str = "pca",
    n_jobs: int = -1,
    n_clusters: int = None,
    random_state: int = 0,
    use_gpu: bool = False,
    max_iter_harmony: int = 10,
) -> str:
    """Batch correction on PCs using Harmony.

    This is a wrapper of `harmony-pytorch <https://github.com/lilab-bcb/harmony-pytorch>`_ package, which is a Pytorch implementation of Harmony algorithm [Korsunsky19]_.

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.

    batch: ``str``, optional, default: ``"Channel"``.
        Which attribute in data.obs field represents batches, default is "Channel".

    rep: ``str``, optional, default: ``"pca"``.
        Which representation to use as input of Harmony, default is PCA.

    n_jobs : ``int``, optional, default: ``-1``.
        Number of threads to use in Harmony. ``-1`` refers to using all physical CPU cores.

    n_clusters: ``int``, optional, default: ``None``.
        Number of Harmony clusters. Default is ``None``, which asks Harmony to estimate this number from the data.

    random_state: ``int``, optional, default: ``0``.
        Seed for random number generator

    use_gpu: ``bool``, optional, default: ``False``.
        If ``True``, use GPU if available. Otherwise, use CPU only.
        
    max_iter_harmony: ``int``, optional, default: ``10``.
        Maximum iterations on running Harmony if not converged.

    Returns
    -------
    out_rep: ``str``
        The keyword in ``data.obsm`` referring to the embedding calculated by Harmony algorithm.

        This keyword is ``rep + '_harmony'``, where ``rep`` is the input parameter above.

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_rep]``: The embedding calculated by Harmony algorithm.

    Examples
    --------
    >>> pg.run_harmony(data, rep = "pca", n_jobs = 10, random_state = 25)
    """
    if not check_batch_key(data, batch, "Cannot apply Harmony!"):
        return rep

    try:
        from harmony import harmonize
    except ImportError as e:
        import sys
        logger.error(f"{e}\nNeed Harmony! Try 'pip install harmony-pytorch'.")
        sys.exit(-1)

    logger.info("Start integration using Harmony.")
    out_rep = rep + "_harmony"
    data.obsm["X_" + out_rep] = harmonize(
        X_from_rep(data, rep),
        data.obs,
        batch,
        n_clusters = n_clusters,
        n_jobs = n_jobs,
        random_state = random_state,
        use_gpu = use_gpu,
        max_iter_harmony = max_iter_harmony,
    )
    return out_rep


@timer(logger=logger)
def run_scanorama(
    data: Union[MultimodalData, UnimodalData],
    batch: str = "Channel",
    n_components: int = 50,
    features: str = "highly_variable_features",
    standardize: bool = True,
    max_value: float = 10.0,
    random_state: int = 0,
) -> str:
    """Batch correction using Scanorama.

    This is a wrapper of `Scanorama <https://github.com/brianhie/scanorama>`_ package. See [Hie19]_ for details on the algorithm.

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.

    batch: ``str``, optional, default: ``"Channel"``.
        Which attribute in data.obs field represents batches, default is "Channel".

    n_components: ``int``, optional default: ``50``.
        Number of integrated embedding components to keep. This sets Scanorama's dimred parameter.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for Scanorama.

    standardize: ``bool``, optional, default: ``True``.
        Whether to scale the data to unit variance and zero mean.

    max_value: ``float``, optional, default: ``10``.
        The threshold to truncate data after scaling. If ``None``, do not truncate.

    random_state: ``int``, optional, default: ``0``.
        Seed for random number generator.

    Returns
    -------
    out_rep: ``str``
        The keyword in ``data.obsm`` referring to the embedding calculated by Scanorama algorithm. out_rep is always equal to "scanorama"

    Update ``data.obsm``:
        * ``data.obsm['X_scanorama']``: The embedding calculated by Scanorama algorithm.

    Examples
    --------
    >>> pg.run_scanorama(data, random_state = 25)
    """
    if not check_batch_key(data, batch, "Cannot apply Scanorama!"):
        return "pca"

    try:
        from scanorama import integrate
    except ImportError as e:
        import sys
        logger.error(f"{e}\nNeed Scanorama! Try 'pip install scanorama'.")
        sys.exit(-1)

    logger.info("Start integration using Scanorama.")

    rep = "scanorama"
    keyword = select_features(data, features=features, standardize=standardize, max_value=max_value)
    X = data.uns[keyword]

    datasets = []
    for channel in data.obs[batch].cat.categories:
        idx = (data.obs[batch] == channel).values
        assert idx.sum() > 0
        datasets.append(X[idx, :])
    genes_list = [[str(i) for i in range(X.shape[1])]] * data.obs[batch].cat.categories.size

    integrated, genes = integrate(datasets, genes_list, dimred = n_components, seed = random_state)
    data.obsm[f"X_{rep}"] = np.concatenate(integrated, axis = 0)

    return rep
