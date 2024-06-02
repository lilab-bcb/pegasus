import time
import numpy as np
import pandas as pd
from anndata import AnnData
from pandas.api.types import is_categorical_dtype
from pegasusio import MultimodalData, UnimodalData
from natsort import natsorted

from threadpoolctl import threadpool_limits
from scipy.sparse import issparse
from sklearn.cluster import KMeans
from typing import List, Optional, Union

from pegasus.tools import eff_n_jobs, construct_graph, calc_stat_per_batch, update_rep, X_from_rep, slicing

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


def _calc_trans_distor(X: np.ndarray, labels: np.ndarray, Y: float) -> float:
    """ Calculate transformed distortion function for the jump method  """
    _, means, _ = calc_stat_per_batch(X, labels)
    distor = (((X - means.T[labels,:]) ** 2).sum() / X.shape[0] / X.shape[1])
    trans_distor = distor ** (-Y)
    return trans_distor

@timer(logger=logger)
def jump_method(
    data: MultimodalData,
    rep: str = "pca",
    K_max: int = 40,
    Y: float = None,
    n_jobs: int = -1,
    random_state: int = 0,
) -> None:
    """ Determine the optimal number of clusters using the Jump Method. [Sugar and James, 2003]_

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

    K_max: ``int``, optional, default: 40
        The maximum number of clusters to try.

    Y: ``float``, optional, default: ``None``
        The transformation power used. If None, use min(data.shape[1] / 3.0, 3.0).

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    Returns
    -------
    ``None``

    Update ``data.uns``:
        * ``data.obs[jump_values]``: Jump values (difference of adjacent transformed distortion values)

    Examples
    --------
    >>> pg.jump_method(data)
    """
    X = data.obsm[f"X_{rep}"]
    Y = min(data.shape[1] / 3.0, 3.0) if Y is None else Y
    logger.info(f"Jump method: Y = {Y:.3f}.")

    n_jobs = eff_n_jobs(n_jobs)
    jump_values = np.zeros(K_max, dtype = np.float64)
    v_old = v = 0.0
    for k in range(1, K_max + 1):
        with threadpool_limits(limits = n_jobs):
            kmeans = KMeans(n_clusters = k, random_state = random_state).fit(X)
        v = _calc_trans_distor(X, kmeans.labels_, Y)
        jump_values[k - 1] = v - v_old
        v_old = v
        logger.info(f"K = {k} is finished, jump_value = {jump_values[k - 1]:.6f}.")
    optimal_k = np.argmax(jump_values) + 1

    data.uns[f"{rep}_jump_values"] = jump_values
    data.uns[f"{rep}_optimal_k"] = optimal_k

    logger.info(f"Jump method finished. Optimal K = {optimal_k}.")


def _run_community_detection(algo, module, G, resolution, random_state, n_iter = -1):
    partition_type = module.RBConfigurationVertexPartition
    if algo == "louvain":
        partition = partition_type(G, resolution_parameter=resolution, weights="weight")
        optimiser = module.Optimiser()
        optimiser.set_rng_seed(random_state)
        diff = optimiser.optimise_partition(partition)
    else:
        partition = module.find_partition(
            G,
            partition_type,
            seed=random_state,
            weights="weight",
            resolution_parameter=resolution,
            n_iterations=n_iter,
        )
    return partition.membership

def _find_optimal_resolution(algo, module, optimal_k, resol_max, G, random_state, n_iter = -1):
    resol = None
    membership = None

    resol_l = 0.01
    resol_r = resol_max
    while resol_r - resol_l > 0.05:
        resol_mid = (resol_l + resol_r) / 2.0
        membership_mid = _run_community_detection(algo, module, G, resol_mid, random_state, n_iter)
        k = max(membership_mid) + 1
        logger.info(f"_find_optimal_resolution: resol = {resol_mid:.4f}, k = {k}, optimal_k = {optimal_k}.")
        if k >= optimal_k:
            resol_r = resol_mid
            resol = resol_mid
            membership = membership_mid
        else:
            resol_l = resol_mid

    if resol is None:
        resol = resol_r
        membership = _run_community_detection(algo, module, G, resol, random_state, n_iter)
        k = max(membership_mid) + 1
        logger.info(f"_find_optimal_resolution: resol = {resol:.4f}, k = {k}, optimal_k = {optimal_k}.")

    return resol, membership


@timer(logger=logger)
def louvain(
    data: MultimodalData,
    rep: str = "pca",
    resolution: int = 1.3,
    n_clust: int = None,
    random_state: int = 0,
    class_label: str = "louvain_labels",
) -> None:
    """Cluster the cells using Louvain algorithm. [Blondel08]_

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm`` and nearest neighbors must be calculated so that affinity matrix ``'W_' + rep`` exists in ``data.uns``. By default, use PCA coordinates.

    resolution: ``int``, optional, default: ``1.3``
        Resolution factor. Higher resolution tends to find more clusters with smaller sizes.

    n_clust: ``int``, optional, default: ``None``
        This option only takes effect if 'resolution = None'. Try to find an appropriate resolution by binary search such that the total number of clusters matches 'n_clust'. The range of resolution to search is (0.01, 2.0].

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    class_label: ``str``, optional, default: ``"louvain_labels"``
        Key name for storing cluster labels in ``data.obs``.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[class_label]``: Cluster labels of cells as categorical data.

    Examples
    --------
    >>> pg.louvain(data)
    """
    try:
        import louvain as louvain_module
    except ImportError:
        import sys
        logger.error("Need louvain! Try 'pip install louvain' or 'conda install -c conda-forge louvain'.")
        sys.exit(-1)
    rep_key = "W_" + rep
    if rep_key not in data.obsp:
        raise ValueError("Cannot find affinity matrix. Please run neighbors first!")
    W = data.obsp[rep_key]

    G = construct_graph(W)
    if resolution is not None:
        membership = _run_community_detection("louvain", louvain_module, G, resolution, random_state)
    else:
        assert isinstance(n_clust, int)
        resolution, membership = _find_optimal_resolution("louvain", louvain_module, n_clust, 2.0, G, random_state)

    data.uns["louvain_resolution"] = resolution
    labels = np.array([str(x + 1) for x in membership])
    categories = natsorted(np.unique(labels))
    data.obs[class_label] = pd.Categorical(values=labels, categories=categories)
    data.register_attr(class_label, "cluster")

    n_clusters = data.obs[class_label].cat.categories.size
    logger.info(f"Louvain clustering is done. Get {n_clusters} clusters.")


@timer(logger=logger)
def leiden(
    data: MultimodalData,
    rep: str = "pca",
    resolution: int = 1.3,
    n_clust: int = None,
    n_iter: int = -1,
    random_state: int = 0,
    class_label: str = "leiden_labels",
) -> None:
    """Cluster the data using Leiden algorithm. [Traag19]_

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm`` and nearest neighbors must be calculated so that affinity matrix ``'W_' + rep`` exists in ``data.uns``. By default, use PCA coordinates.

    resolution: ``int``, optional, default: ``1.3``
        Resolution factor. Higher resolution tends to find more clusters.

    n_clust: ``int``, optional, default: ``None``
        This option only takes effect if 'resolution = None'. Try to find an appropriate resolution by binary search such that the total number of clusters matches 'n_clust'. The range of resolution to search is (0.01, 2.0].

    n_iter: ``int``, optional, default: ``-1``
        Number of iterations that Leiden algorithm runs. If ``-1``, run the algorithm until reaching its optimal clustering.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    class_label: ``str``, optional, default: ``"leiden_labels"``
        Key name for storing cluster labels in ``data.obs``.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[class_label]``: Cluster labels of cells as categorical data.

    Examples
    --------
    >>> pg.leiden(data)
    """

    try:
        import leidenalg
    except ImportError:
        import sys
        logger.error("Need leidenalg! Try 'pip install leidenalg'.")
        sys.exit(-1)

    rep_key = "W_" + rep
    if rep_key not in data.obsp:
        raise ValueError("Cannot find affinity matrix. Please run neighbors first!")
    W = data.obsp[rep_key]

    G = construct_graph(W)
    if resolution is not None:
        membership = _run_community_detection("leiden", leidenalg, G, resolution, random_state, n_iter)
    else:
        assert isinstance(n_clust, int)
        resolution, membership = _find_optimal_resolution("leiden", leidenalg, n_clust, 2.0, G, random_state, n_iter)

    data.uns["leiden_resolution"] = resolution
    labels = np.array([str(x + 1) for x in membership])
    categories = natsorted(np.unique(labels))
    data.obs[class_label] = pd.Categorical(values=labels, categories=categories)
    data.register_attr(class_label, "cluster")

    n_clusters = data.obs[class_label].cat.categories.size
    logger.info(f"Leiden clustering is done. Get {n_clusters} clusters.")


def partition_cells_by_kmeans(
    X: np.ndarray,
    n_clusters: int,
    n_clusters2: int,
    n_init: int,
    n_jobs: int,
    random_state: int,
    min_avg_cells_per_final_cluster: Optional[int] = 10,
) -> List[int]:

    n_clusters = min(n_clusters, max(X.shape[0] // min_avg_cells_per_final_cluster, 1))
    if n_clusters == 1:
        return np.zeros(X.shape[0], dtype = np.int32)

    n_jobs = eff_n_jobs(n_jobs)

    kmeans_params = {
        'n_clusters': n_clusters,
        'n_init': n_init,
        'random_state': random_state,
    }
    km = KMeans(**kmeans_params)

    with threadpool_limits(limits = n_jobs):
        km.fit(X)
        coarse = km.labels_.copy()

        km.set_params(n_init=1)
        labels = coarse.copy()
        base_sum = 0
        for i in range(n_clusters):
            idx = coarse == i
            nc = min(n_clusters2, max(idx.sum() // min_avg_cells_per_final_cluster, 1))
            if nc == 1:
                labels[idx] = base_sum
            else:
                km.set_params(n_clusters=nc)
                km.fit(X[idx, :])
                labels[idx] = base_sum + km.labels_
            base_sum += nc

    return labels


@timer(logger=logger)
def spectral_louvain(
    data: MultimodalData,
    rep: str = "pca",
    resolution: float = 1.3,
    rep_kmeans: str = "diffmap",
    n_clusters: int = 30,
    n_clusters2: int = 50,
    n_init: int = 10,
    n_jobs: int = -1,
    random_state: int = 0,
    class_label: str = "spectral_louvain_labels",
) -> None:
    """ Cluster the data using Spectral Louvain algorithm. [Li20]_

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

    resolution: ``int``, optional, default: ``1.3``
        Resolution factor. Higher resolution tends to find more clusters with smaller sizes.

    rep_kmeans: ``str``, optional, default: ``"diffmap"``
        The embedding representation on which the KMeans runs. Keyword must exist in ``data.obsm``. By default, use Diffusion Map coordinates. If diffmap is not calculated, use PCA coordinates instead.

    n_clusters: ``int``, optional, default: ``30``
        The number of first level clusters.

    n_clusters2: ``int``, optional, default: ``50``
        The number of second level clusters.

    n_init: ``int``, optional, default: ``10``
        Number of kmeans tries for the first level clustering. Default is set to be the same as scikit-learn Kmeans function.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use for the KMeans step. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    class_label: ``str``, optional, default: ``"spectral_louvain_labels"``
        Key name for storing cluster labels in ``data.obs``.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[class_label]``: Cluster labels for cells as categorical data.

    Examples
    --------
    >>> pg.spectral_louvain(data)
    """
    try:
        import louvain as louvain_module
    except ImportError:
        import sys
        logger.error("Need louvain! Try 'pip install louvain' or 'conda install -c conda-forge louvain'.")
        sys.exit(-1)

    if f"X_{rep_kmeans}" not in data.obsm.keys():
        logger.warning(f"{rep_kmeans} is not calculated, switch to pca instead.")
        rep_kmeans = "pca"
        if f"X_{rep_kmeans}" not in data.obsm.keys():
            raise ValueError(f"Please run {rep_kmeans} first!")
    if f"W_{rep}" not in data.obsp:
        raise ValueError("Cannot find affinity matrix. Please run neighbors first!")

    labels = partition_cells_by_kmeans(
        data.obsm[f"X_{rep_kmeans}"], n_clusters, n_clusters2, n_init, n_jobs, random_state,
    )

    W = data.obsp[f"W_{rep}"]

    G = construct_graph(W)
    partition_type = louvain_module.RBConfigurationVertexPartition
    partition = partition_type(
        G, resolution_parameter=resolution, weights="weight", initial_membership=labels
    )
    partition_agg = partition.aggregate_partition()

    optimiser = louvain_module.Optimiser()
    optimiser.set_rng_seed(random_state)
    diff = optimiser.optimise_partition(partition_agg)
    partition.from_coarse_partition(partition_agg)

    labels = np.array([str(x + 1) for x in partition.membership])
    categories = natsorted(np.unique(labels))
    data.obs[class_label] = pd.Categorical(values=labels, categories=categories)
    data.register_attr(class_label, "cluster")

    n_clusters = data.obs[class_label].cat.categories.size
    logger.info(f"Spectral Louvain clustering is done. Get {n_clusters} clusters.")


@timer(logger=logger)
def spectral_leiden(
    data: MultimodalData,
    rep: str = "pca",
    resolution: float = 1.3,
    rep_kmeans: str = "diffmap",
    n_clusters: int = 30,
    n_clusters2: int = 50,
    n_init: int = 10,
    n_jobs: int = -1,
    random_state: int = 0,
    class_label: str = "spectral_leiden_labels",
) -> None:
    """Cluster the data using Spectral Leiden algorithm. [Li20]_

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

    resolution: ``int``, optional, default: ``1.3``
        Resolution factor. Higher resolution tends to find more clusters.

    rep_kmeans: ``str``, optional, default: ``"diffmap"``
        The embedding representation on which the KMeans runs. Keyword must exist in ``data.obsm``. By default, use Diffusion Map coordinates. If diffmap is not calculated, use PCA coordinates instead.

    n_clusters: ``int``, optional, default: ``30``
        The number of first level clusters.

    n_clusters2: ``int``, optional, default: ``50``
        The number of second level clusters.

    n_init: ``int``, optional, default: ``10``
        Number of kmeans tries for the first level clustering. Default is set to be the same as scikit-learn Kmeans function.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use for the KMeans step. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    class_label: ``str``, optional, default: ``"spectral_leiden_labels"``
        Key name for storing cluster labels in ``data.obs``.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[class_label]``: Cluster labels for cells as categorical data.

    Examples
    --------
    >>> pg.spectral_leiden(data)
    """

    try:
        import leidenalg
    except ImportError:
        import sys
        logger.error("Need leidenalg! Try 'pip install leidenalg'.")
        sys.exit(-1)

    if f"X_{rep_kmeans}" not in data.obsm.keys():
        logger.warning(f"{rep_kmeans} is not calculated, switch to pca instead.")
        rep_kmeans = "pca"
        if f"X_{rep_kmeans}" not in data.obsm.keys():
            raise ValueError(f"Please run {rep_kmeans} first!")
    if f"W_{rep}" not in data.obsp:
        raise ValueError("Cannot find affinity matrix. Please run neighbors first!")

    labels = partition_cells_by_kmeans(
        data.obsm[f"X_{rep_kmeans}"], n_clusters, n_clusters2, n_init, n_jobs, random_state,
    )

    W = data.obsp[f"W_{rep}"]

    G = construct_graph(W)
    partition_type = leidenalg.RBConfigurationVertexPartition
    partition = partition_type(
        G, resolution_parameter=resolution, weights="weight", initial_membership=labels
    )
    partition_agg = partition.aggregate_partition()

    optimiser = leidenalg.Optimiser()
    optimiser.set_rng_seed(random_state)
    diff = optimiser.optimise_partition(partition_agg, -1)
    partition.from_coarse_partition(partition_agg)

    labels = np.array([str(x + 1) for x in partition.membership])
    categories = natsorted(np.unique(labels))
    data.obs[class_label] = pd.Categorical(values=labels, categories=categories)
    data.register_attr(class_label, "cluster")

    n_clusters = data.obs[class_label].cat.categories.size
    logger.info(f"Spectral Leiden clustering is done. Get {n_clusters} clusters.")


def cluster(
    data: MultimodalData,
    algo: str = "louvain",
    rep: str = "pca",
    resolution: int = 1.3,
    n_jobs: int = -1,
    random_state: int = 0,
    class_label: str = None,
    n_iter: int = -1,
    rep_kmeans: str = "diffmap",
    n_clusters: int = 30,
    n_clusters2: int = 50,
    n_init: int = 10,
) -> None:
    """Cluster the data using the chosen algorithm.

    Candidates are *louvain*, *leiden*, *spectral_louvain* and *spectral_leiden*.
    If data have < 1000 cells and there are clusters with sizes of 1, resolution is automatically reduced until no cluster of size 1 appears.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    algo: ``str``, optional, default: ``"louvain"``
        Which clustering algorithm to use. Choices are louvain, leiden, spectral_louvain, spectral_leiden

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

    resolution: ``int``, optional, default: ``1.3``
        Resolution factor. Higher resolution tends to find more clusters.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use for the KMeans step in 'spectral_louvain' and 'spectral_leiden'. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    class_label: ``str``, optional, default: None
        Key name for storing cluster labels in ``data.obs``. If None, use 'algo_labels'.

    n_iter: ``int``, optional, default: ``-1``
        Number of iterations that Leiden algorithm runs. If ``-1``, run the algorithm until reaching its optimal clustering.

    rep_kmeans: ``str``, optional, default: ``"diffmap"``
        The embedding representation on which the KMeans runs. Keyword must exist in ``data.obsm``. By default, use Diffusion Map coordinates. If diffmap is not calculated, use PCA coordinates instead.

    n_clusters: ``int``, optional, default: ``30``
        The number of first level clusters.

    n_clusters2: ``int``, optional, default: ``50``
        The number of second level clusters.

    n_init: ``int``, optional, default: ``10``
        Number of kmeans tries for the first level clustering. Default is set to be the same as scikit-learn Kmeans function.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[class_label]``: Cluster labels of cells as categorical data.

    Examples
    --------
    >>> pg.cluster(data, algo = 'leiden')
    """

    if algo not in {"louvain", "leiden", "spectral_louvain", "spectral_leiden"}:
        raise ValueError("Unknown clustering algorithm {}.".format(algo))

    if class_label is None:
        class_label = algo + "_labels"

    kwargs = {
        "data": data,
        "rep": rep,
        "resolution": resolution,
        "random_state": random_state,
        "class_label": class_label,
    }
    if algo == "leiden":
        kwargs["n_iter"] = n_iter
    if algo in ["spectral_louvain", "spectral_leiden"]:
        kwargs.update(
            {
                "rep_kmeans": rep_kmeans,
                "n_clusters": n_clusters,
                "n_clusters2": n_clusters2,
                "n_init": n_init,
                "n_jobs": n_jobs,
            }
        )

    cluster_func = globals()[algo]

    cluster_func(**kwargs)  # clustering
    if data.shape[0] < 100000 and data.obs[class_label].value_counts().min() == 1:
        new_resol = resolution
        while new_resol > 0.0:
            new_resol -= 0.1
            kwargs["resolution"] = new_resol
            cluster_func(**kwargs)
            if data.obs[class_label].value_counts().min() > 1:
                break
        logger.warning(
            "Reduced resolution from {:.2f} to {:.2f} to avoid clusters of size 1.".format(
                resolution, new_resol
            )
        )


def split_one_cluster(
    data: MultimodalData,
    clust_label: str,
    clust_id: str,
    n_clust: int,
    res_label: str,
    rep: str = "pca",
    n_comps: int = None,
    random_state: int = 0,
) -> None:
    """
    Use Leiden algorithm to split 'clust_id' in 'clust_label' into 'n_components' sub-clusters and write the new clusting results to 'res_label'. The sub-cluster names are the concatenation of original cluster name and the subcluster id (e.g. 'T' -> 'T-1', 'T-2').

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    clust_label: `str`
        Use existing clustering stored in data.obs['clust_label'].

    clust_id: `str`
        Cluster ID in data.obs['clust_label'].

    n_clust: `int`
        Split 'clust_id' into `n_clust' subclusters.

    res_label: `str`,
        Write new clustering in data.obs['res_label']. The sub-cluster names are the concatenation of original cluster name and the subcluster id (e.g. 'T' -> 'T-1', 'T-2').

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for Kmeans clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

    n_comps: `int`, optional (default: None)
        Number of components to be used in the `rep`. If n_comps == None, use all components; otherwise, use the minimum of n_comps and rep's dimensions.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use for the KMeans step in 'spectral_louvain' and 'spectral_leiden'. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[res_label]``: New cluster labels of cells as categorical data.

    Examples
    --------
    >>> pg.split_one_cluster(data, 'leiden_labels', '15', 2, 'leiden_labels_split')
    """
    cats = None
    if is_categorical_dtype(data.obs[clust_label]):
        cats = data.obs[clust_label].cat.categories.values
    else:
        cats = pd.Categorical(data.obs[clust_label]).categories.values
        if cats.dtype.kind not in {'S', 'U'}:
            cats = cats.astype(str)
    idx_cat = np.nonzero(cats==clust_id)[0]

    if idx_cat.size == 0:
        raise ValueError(f"{clust_id} is not in {clust_label}!")
    elif idx_cat.size > 1:
        raise ValueError(f"Detected more than one categories in {clust_label} with name {clust_id}!")
    else:
        idx_cat = idx_cat[0]

    idx = np.nonzero((data.obs[clust_label] == clust_id).values)[0]
    tmpdat = data[idx].copy()
    from pegasus.tools import neighbors
    neighbors(tmpdat, rep=rep, n_comps=n_comps, use_cache=False)
    leiden(tmpdat, rep=rep, resolution=None, n_clust=n_clust, random_state=random_state)

    new_clust = data.obs[clust_label].values.astype(object)
    cats_sub = []
    for i, label in enumerate(tmpdat.obs['leiden_labels'].value_counts().index):
        sub_id = f"{clust_id}-{i+1}"
        new_clust[idx[(tmpdat.obs['leiden_labels'] == label).values]] = sub_id
        cats_sub.append(sub_id)

    data.obs[res_label] = pd.Categorical(values = new_clust, categories = np.concatenate((cats[0:idx_cat], np.array(cats_sub), cats[idx_cat+1:])))
    data.register_attr(res_label, "cluster")
    del tmpdat


@timer(logger=logger)
def calc_dendrogram(
    data: Union[MultimodalData, UnimodalData, AnnData],
    groupby: str = "obs",
    rep: Optional[str] = "pca",
    genes: Optional[List[str]] = None,
    on_average: bool = True,
    linkage_method: str = "ward",
    res_key: str = "dendrogram",
) -> None:
    """
    Cluster data using hierarchical clustering algorithm.

    The metric in use is a Connection Specific Index (CSI) matrix ([Suo18]_, [Bass13]_) built from the correlations between ``groupby`` attribute levels regarding the ``rep`` embedding.

    Parameters
    ----------

    data: ``MultimodalData``, ``UnimodalData``, or ``AnnData`` object
        Single cell expression data.
    groupby: ``str``, optional, default: ``None``
        Set cluster labels in use.
        If ``"obs"``, use cell names (i.e. ``data.obs_names``); if ``"var"``, use feature names (i.e. ``data.var_names``).
        Otherwise, specify a categorical cell or feature attribute to use, which must exist in ``data.obs`` or ``data.var``.
    rep: ``str``, optional, default: ``pca``
        Cell embedding to use. If specified, it only works when ``genes`` is ``None``, and its key ``"X_"+rep`` must exist in ``data.obsm``. By default, use PCA embedding.
        If ``None``, use the current count matrix ``data.X``.
    genes: ``List[str]``, optional, default: ``None``
        List of genes to use. Gene names must exist in ``data.var``. If set, use the counts in ``data.X`` for plotting; if ``None``, use the embedding specified in ``rep``.
    on_average: ``bool``, optional, default: ``True``
        If ``True``, clustering ``groupby`` levels based on their mean values. Only works when ``groupby`` is not ``None``.
    linkage_method: ``str``, optional, default: ``ward``
        Which linkage criterion to use, used by hierarchical clustering. Available options: ``ward`` (default), ``single``, ``complete``, ``average``, ``weighted``, ``centroid``, ``median``.
        See `scipy linkage documentation`_ for details.
    res_key: ``str``, optional, default: ``dendrogram``
        Key name in ``data.uns`` field to store the calculated dendrogram information, which will be used by ``plot_dendrogram`` function for plotting.

    Returns
    -------
    ``None``

    Update ``data.uns``:
        * ``data.uns[res_key]``: A tuple of the calculated linkage matrix and its corresponding labels.

    Examples
    --------
    >>> pg.calc_dendrogram(data, groupby='leiden_labels')
    >>> pg.calc_dendrogram(data, genes=['CD4', 'CD8A', 'CD8B'], on_average=False)
    >>> pg.calc_dendrogram(data, groupby="var", rep=None, on_average=False)

    .. _scipy linkage documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    """
    # Set up embedding or count matrix to use
    if genes is None:
        rep = update_rep(rep)
        embed_df = pd.DataFrame(X_from_rep(data, rep))
    else:
        embed_df = pd.DataFrame(slicing(data[:, genes].X))

    # Set up index
    if groupby == "obs":
        indices = data.obs_names
    elif groupby == "var":
        embed_df = embed_df.T
        indices = data.var_names
    elif groupby in data.obs:
        indices = data.obs[groupby]
    elif groupby in data.var:
        embed_df = embed_df.T
        indices = data.var[groupby]
    else:
        raise Exception(f"The groupby key {groupby} doesn't exist in data.obs or data.var!")
    embed_df.set_index(indices, inplace=True)

    # Use group mean if on_average is True
    if on_average:
        embed_df = embed_df.groupby(level=0, observed=True).mean()
    if not isinstance(embed_df.index.dtype, pd.CategoricalDtype):
        embed_df.index = embed_df.index.astype("category")

    # Calculate Pearson's correlation between cluster labels
    corr_df = pd.DataFrame(np.corrcoef(embed_df, rowvar=True), columns=embed_df.index, index=embed_df.index)  # Faster than pandas corr
    corr_mat = corr_df.values

    from pegasus.tools.utils import calc_csi_matrix

    # Calculate CSI matrix
    csi_mat = calc_csi_matrix(corr_mat)
    csi_df = pd.DataFrame(csi_mat, columns=corr_df.index, index=corr_df.index)

    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform

    dissim_df = 1 - csi_df
    np.fill_diagonal(dissim_df.to_numpy(), 0)    # Enforce main diagonal to be 0 to pass squareform requirement
    Z = linkage(squareform(dissim_df), method=linkage_method, optimal_ordering=True)

    data.uns[res_key] = (Z, csi_df)
