import time
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from scipy.sparse import issparse
from pegasusio import MultimodalData

from pegasus.tools import estimate_feature_statistics, select_features, X_from_rep

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



def set_group_attribute(data: MultimodalData, attribute_string: str) -> None:
    """Set group attributes used in batch correction.

    Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, *pegasus* will only perform batch correction for channels within each group.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    attribute_string: ``str``
        Attributes used to construct groups:

        * ``None``
            Assume all channels are from one group.

        * ``attr``
            Define groups by sample attribute ``attr``, which is a keyword in ``data.obs``.

        * ``att1+att2+...+attrn``
            Define groups by the Cartesian product of these *n* attributes, which are keywords in ``data.obs``.

        * ``attr=value_11,...value_1n_1;value_21,...value_2n_2;...;value_m1,...,value_mn_m``
            In this form, there will be *(m+1)* groups. A cell belongs to group *i* (*i > 1*) if and only if its sample attribute ``attr``, which is a keyword in ``data.obs``, has a value among ``value_i1``, ... ``value_in_i``. A cell belongs to group 0 if it does not belong to any other groups.

    Returns
    -------
    None

        Update ``data.obs``:

        * ``data.obs["Group"]``: Group ID for each cell.

    Examples
    --------

    >>> pg.set_group_attribute(data, attr_string = "Individual")

    >>> pg.set_group_attribute(data, attr_string = "Individual+assignment")

    >>> pg.set_group_attribute(data, attr_string = "Channel=1,3,5;2,4,6,8")
    """

    if attribute_string.find("=") >= 0:
        attr, value_str = attribute_string.split("=")
        assert attr in data.obs.columns
        values = value_str.split(";")
        data.obs["Group"] = "0"
        for group_id, value in enumerate(values):
            vals = value.split(",")
            idx = np.isin(data.obs[attr], vals)
            data.obs.loc[idx, "Group"] = str(group_id + 1)
    elif attribute_string.find("+") >= 0:
        attrs = attribute_string.split("+")
        assert np.isin(attrs, data.obs.columns).sum() == len(attrs)
        data.obs["Group"] = data.obs[attrs].apply(lambda x: "+".join(x), axis=1)
    else:
        assert attribute_string in data.obs.columns
        data.obs["Group"] = data.obs[attribute_string]


def estimate_adjustment_matrices(data: MultimodalData) -> bool:
    """ Estimate adjustment matrices
    """

    if "plus" in data.varm.keys() or "muls" in data.varm.keys():
        # This only happens if this is for subclustering. Thus do not calculate factors, using factors calculated from parent for batch correction.
        assert "plus" in data.varm.keys() and "muls" in data.varm.keys()
        return True

    if ("gmeans" not in data.varm) or ("gstds" not in data.varm):
        estimate_feature_statistics(data, True)

    if data.uns["Channels"].size == 1:
        logger.warning(
            "Warning: data only contains 1 channel. Batch correction disabled!"
        )
        return False

    nchannel = data.uns["Channels"].size

    plus = np.zeros((data.shape[1], nchannel))
    muls = np.zeros((data.shape[1], nchannel))

    ncells = data.uns["ncells"]
    means = data.varm["means"]
    partial_sum = data.varm["partial_sum"]
    gmeans = data.varm["gmeans"]
    gstds = data.varm["gstds"]
    c2gid = data.uns["c2gid"]
    for i in range(data.uns["Channels"].size):
        if ncells[i] > 1:
            muls[:, i] = (partial_sum[:, i] / (ncells[i] - 1.0)) ** 0.5
        outliers = muls[:, i] < 1e-6
        normals = np.logical_not(outliers)
        muls[outliers, i] = 1.0
        muls[normals, i] = gstds[normals, c2gid[i]] / muls[normals, i]
        plus[:, i] = gmeans[:, c2gid[i]] - muls[:, i] * means[:, i]

    data.varm["plus"] = plus
    data.varm["muls"] = muls

    return True

def correct_batch_effects(data: MultimodalData, keyword: str, features: str = None) -> None:
    """ Apply calculated plus and muls to correct batch effects for a dense matrix
    """
    X = data.uns[keyword]
    m = X.shape[1]
    if features is not None:
        selected = data.var[features].values
        plus = data.varm["plus"][selected, :]
        muls = data.varm["muls"][selected, :]
    else:
        selected = np.ones(data.shape[1], dtype=bool)
        plus = data.varm["plus"]
        muls = data.varm["muls"]

    for i, channel in enumerate(data.uns["Channels"]):
        idx = np.isin(data.obs["Channel"], channel)
        if idx.sum() == 0:
            continue
        X[idx] = X[idx] * np.reshape(muls[:, i], newshape=(1, m)) + np.reshape(
            plus[:, i], newshape=(1, m)
        )
    # X[X < 0.0] = 0.0


def correct_batch(data: MultimodalData, features: str = None) -> None:
    """Batch correction on data.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    features: `str`, optional, default: ``None``
        Features to be included in batch correction computation. If ``None``, simply consider all features.

    Returns
    -------
    ``None``

    Update ``data.X`` by the corrected count matrix.

    Examples
    --------
    >>> pg.correct_batch(data, features = "highly_variable_features")
    """

    tot_seconds = 0.0

    # estimate adjustment parameters
    start = time.perf_counter()
    can_correct = estimate_adjustment_matrices(data)
    end = time.perf_counter()
    tot_seconds += end - start
    logger.info("Adjustment parameters are estimated.")

    # select dense matrix
    keyword = select_features(data, features = features, standardize = False, max_value = None) # do not standardize or truncate max_value
    logger.info("Features are selected.")

    if can_correct:
        start = time.perf_counter()
        correct_batch_effects(data, keyword, features)
        end = time.perf_counter()
        tot_seconds += end - start
        logger.info(
            "Batch correction is finished. Time spent = {:.2f}s.".format(tot_seconds)
        )


@timer(logger=logger)
def run_harmony(
    data: MultimodalData,
    rep: str = 'pca',
    n_jobs: int = -1,
    n_clusters: int = None,
    random_state: int = 0,
) -> str:
    """Batch correction PCs using Harmony

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``.
        Which representation to use as input of Harmony, default is PCA.

    n_jobs : ``int``, optional, default: ``-1``.
        Number of threads to use for the KMeans clustering used in Harmony. ``-1`` refers to using all available threads.

    n_clusters: ``int``, optional, default: ``None``.
        Number of Harmony clusters. Default is ``None``, which asks Harmony to estimate this number from the data.

    random_state: ``int``, optional, default: ``0``.
        Seed for random number generator

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
    if not is_categorical_dtype(data.obs['Channel']):
        data.obs['Channel'] = pd.Categorical(data.obs['Channel'])
    if data.obs['Channel'].cat.categories.size  == 1:
        logger.warning("Warning: data only contains 1 channel. Cannot apply Harmony!")
        return rep

    try:
        from harmony import harmonize
    except ImportError as e:
        print(f"ERROR: {e}")
        print("ERROR: Need Harmony! Try 'pip install harmony-pytorch'.")
        import sys
        sys.exit(-1)


    logger.info("Start integration using Harmony.")
    out_rep = rep + '_harmony'
    data.obsm['X_' + out_rep] = harmonize(X_from_rep(data, rep), data.obs, 'Channel', n_clusters = n_clusters, n_jobs_kmeans = n_jobs, random_state = random_state)
    return out_rep


@timer(logger=logger)
def run_scanorama(
    data: MultimodalData,
    n_components: int = 50,
    features: str = "highly_variable_features",
    standardize: bool = True,
    max_value: float = 10,
    random_state: int = 0,
) -> str:
    """Batch correction using Scanorama

    Parameters
    ----------
    data: ``MultimodalData``.
        Annotated data matrix with rows for cells and columns for genes.
    
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
    if not is_categorical_dtype(data.obs['Channel']):
        data.obs['Channel'] = pd.Categorical(data.obs['Channel'])
    if data.obs['Channel'].cat.categories.size  == 1:
        logger.warning("Warning: data only contains 1 channel. Cannot apply Scanorama!")
        return 'pca'

    try:
        from scanorama import integrate
    except ImportError as e:
        print(f"ERROR: {e}")
        print("ERROR: Need Scanorama! Try 'pip install scanorama'.")
        import sys
        sys.exit(-1)

    logger.info("Start integration using Scanorama.")

    rep = 'scanorama'
    keyword = select_features(data, features = features, standardize = standardize, max_value = max_value)
    X = data.uns[keyword]

    datasets = []
    for channel in data.obs['Channel'].cat.categories:
        idx = (data.obs['Channel'] == channel).values
        assert idx.sum() > 0
        datasets.append(X[idx, :])
    genes_list = [[str(i) for i in range(X.shape[1])]] * data.obs['Channel'].cat.categories.size
    
    integrated, genes = integrate(datasets, genes_list, dimred = n_components, seed = random_state)
    data.obsm[f'X_{rep}'] = np.concatenate(integrated, axis = 0)

    return rep
