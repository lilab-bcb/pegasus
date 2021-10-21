import time
import numpy as np
import pandas as pd
from pegasusio import MultimodalData
from natsort import natsorted

from threadpoolctl import threadpool_limits
from sklearn.cluster import KMeans
from typing import List, Optional, Union

from pegasus.tools import eff_n_jobs, construct_graph, calc_stat_per_batch

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
        print("Need louvain! Try 'pip install louvain-github'.")

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
        print("Need leidenalg! Try 'pip install leidenalg'.")

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
        print("Need louvain! Try 'pip install louvain-github'.")

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
        print("Need leidenalg! Try 'pip install leidenalg'.")

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
    random_state: int = 0,
) -> None:
    """
    Use Leiden algorithm to split 'clust_id' in 'clust_label' into 'n_components' clusters and write the new clusting results to 'res_label'. Assume 'clust_label' named clusters as numbers (in str format).

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
        Write new clustering in data.obs['res_label']. The largest subcluster will use 'clust_id' as its cluster ID, while other subclusters will be numbered after existing clusters.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation used for Kmeans clustering. Keyword ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

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
    idx = np.where(data.obs[clust_label] == clust_id)[0]
    tmpdat = data[idx].copy()
    from pegasus.tools import neighbors
    neighbors(tmpdat, rep=rep)
    leiden(tmpdat, rep=rep, resolution=None, n_clust=n_clust, random_state=random_state)
    new_clust = data.obs[clust_label].values.astype(int)
    new_label = new_clust.max() + 1
    for label in tmpdat.obs['leiden_labels'].value_counts().index[1:]:
        new_clust[idx[(tmpdat.obs['leiden_labels'] == label).values]] = new_label
        new_label += 1
    data.obs[res_label] = pd.Categorical(values = new_clust.astype(str), categories = np.array(range(1, new_label)).astype(str))
    data.register_attr(res_label, "cluster")
    del tmpdat
