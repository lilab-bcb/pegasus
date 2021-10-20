import time
import numpy as np
import pandas as pd
import numba

from scipy.sparse import issparse, csr_matrix
from scipy.stats import chi2
from sklearn.neighbors import NearestNeighbors
from pegasusio import MultimodalData
from typing import List, Tuple

from pegasus.tools import eff_n_jobs, update_rep, X_from_rep

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



@numba.njit(cache=True)
def _reorg_knn(indices, distances):
    for i in range(indices.shape[0]):
        if indices[i, 0] != i:
            for j in range(indices.shape[1]-1, 0, -1):
                indices[i, j] = indices[i, j-1]
                distances[i, j] = distances[i, j-1]
            indices[i, 0] = i
            distances[i, 0] = 0.0


def calculate_nearest_neighbors(
    X: np.array,
    K: int = 100,
    n_jobs: int = -1,
    method: str = "hnsw",
    M: int = 20,
    efC: int = 200,
    efS: int = 200,
    random_state: int = 0,
    full_speed: int = False,
):
    """Calculate nearest neighbors
    X is the sample by feature matrix
    Return K -1 neighbors, the first one is the point itself and thus omitted.
    TODO: Documentation
    """

    nsample = X.shape[0]

    if nsample <= 1000:
        method = "sklearn"

    if nsample < K:
        logger.warning(f"Warning: in calculate_nearest_neighbors, number of samples = {nsample} < K = {K}!\n Set K to {nsample}.")
        K = nsample

    n_jobs = eff_n_jobs(n_jobs)

    if method == "hnsw":
        import hnswlib

        assert not issparse(X)
        # Build hnsw index
        knn_index = hnswlib.Index(space="l2", dim=X.shape[1])
        knn_index.init_index(
            max_elements=nsample, ef_construction=efC, M=M, random_seed=random_state
        )
        knn_index.set_num_threads(n_jobs if full_speed else 1)
        knn_index.add_items(X)

        # KNN query
        knn_index.set_ef(efS)
        knn_index.set_num_threads(n_jobs)
        indices, distances = knn_index.knn_query(X, k=K)
        # eliminate the first neighbor, which is the node itself
        _reorg_knn(indices, distances)
        indices = indices[:, 1:]
        indices.dtype = np.int64
        distances = distances[:, 1:]
        distances = np.sqrt(distances, out=distances)
    else:
        assert method == "sklearn"
        knn = NearestNeighbors(
            n_neighbors=K - 1, n_jobs=n_jobs
        )  # eliminate the first neighbor, which is the node itself
        knn.fit(X)
        distances, indices = knn.kneighbors()

    return indices, distances


def knn_is_cached(
    data: MultimodalData, indices_key: str, distances_key: str, K: int
) -> bool:
    return (
        (indices_key in data.obsm)
        and (distances_key in data.obsm)
        and data.obsm[indices_key].shape[0] == data.shape[0]
        and (K <= data.obsm[indices_key].shape[1] + 1)
    )


@timer(logger=logger)
def get_neighbors(
    data: MultimodalData,
    K: int = 100,
    rep: str = "pca",
    n_jobs: int = -1,
    random_state: int = 0,
    full_speed: bool = False,
    use_cache: bool = True,
) -> Tuple[List[int], List[float]]:
    """Find K nearest neighbors for each data point and return the indices and distances arrays.

    Parameters
    ----------

    data : `pegasusio.MultimodalData`
        An AnnData object.
    K : `int`, optional (default: 100)
        Number of neighbors, including the data point itself.
    rep : `str`, optional (default: 'pca')
        Representation used to calculate kNN. If `None` use data.X
    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.
    random_state: `int`, optional (default: 0)
        Random seed for random number generator.
    full_speed: `bool`, optional (default: False)
        If full_speed, use multiple threads in constructing hnsw index. However, the kNN results are not reproducible. If not full_speed, use only one thread to make sure results are reproducible.
    use_cache: `bool`, optional (default: True)
        If use_cache and found cached knn results, will not recompute.

    Returns
    -------

    kNN indices and distances arrays.

    Examples
    --------
    >>> indices, distances = tools.get_neighbors(data)
    """

    rep = update_rep(rep)
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if use_cache and knn_is_cached(data, indices_key, distances_key, K):
        indices = data.obsm[indices_key]
        distances = data.obsm[distances_key]
        logger.info("Found cached kNN results, no calculation is required.")
    else:
        indices, distances = calculate_nearest_neighbors(
            X_from_rep(data, rep),
            K=K,
            n_jobs=eff_n_jobs(n_jobs),
            random_state=random_state,
            full_speed=full_speed,
        )
        data.obsm[indices_key] = indices
        data.register_attr(indices_key, "knn")
        data.obsm[distances_key] = distances
        data.register_attr(distances_key, "knn")

    return indices, distances


def get_symmetric_matrix(csr_mat: "csr_matrix") -> "csr_matrix":
    tp_mat = csr_mat.transpose().tocsr()
    sym_mat = csr_mat + tp_mat
    sym_mat.sort_indices()

    idx_mat = (csr_mat != 0).astype(int) + (tp_mat != 0).astype(int)
    idx_mat.sort_indices()
    idx = idx_mat.data == 2

    sym_mat.data[idx] /= 2.0
    return sym_mat


# We should not modify distances array!
@timer(logger=logger)
def calculate_affinity_matrix(
    indices: List[int], distances: List[float]
) -> "csr_matrix":

    nsample = indices.shape[0]
    K = indices.shape[1]
    # calculate sigma, important to use median here!
    sigmas = np.median(distances, axis=1)
    sigmas_sq = np.square(sigmas)

    # calculate local-scaled kernel
    normed_dist = np.zeros((nsample, K), dtype=np.float32)
    for i in range(nsample):
        numers = 2.0 * sigmas[i] * sigmas[indices[i, :]]
        denoms = sigmas_sq[i] + sigmas_sq[indices[i, :]]
        normed_dist[i, :] = np.sqrt(numers / denoms) * np.exp(
            -np.square(distances[i, :]) / denoms
        )

    W = csr_matrix(
        (normed_dist.ravel(), (np.repeat(range(nsample), K), indices.ravel())),
        shape=(nsample, nsample),
    )
    W = get_symmetric_matrix(W)

    # density normalization
    z = W.sum(axis=1).A1
    W = W.tocoo()
    W.data /= z[W.row]
    W.data /= z[W.col]
    W = W.tocsr()
    W.eliminate_zeros()

    return W


def neighbors(
    data: MultimodalData,
    K: int = 100,
    rep: "str" = "pca",
    n_jobs: int = -1,
    random_state: int = 0,
    full_speed: bool = False,
    use_cache: bool = True,
) -> None:
    """Compute k nearest neighbors and affinity matrix, which will be used for diffmap and graph-based community detection algorithms.

    The kNN calculation uses `hnswlib <https://github.com/nmslib/hnswlib>`_ introduced by [Malkov16]_.

    Parameters
    ----------

    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    K: ``int``, optional, default: ``100``
        Number of neighbors, including the data point itself.

    rep: ``str``, optional, default: ``"pca"``
        Embedding representation used to calculate kNN. If ``None``, use ``data.X``; otherwise, keyword ``'X_' + rep`` must exist in ``data.obsm``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    full_speed: ``bool``, optional, default: ``False``
        * If ``True``, use multiple threads in constructing ``hnsw`` index. However, the kNN results are not reproducible.
        * Otherwise, use only one thread to make sure results are reproducible.

    use_cache: ``bool``, optional, default: ``True``
        * If ``True`` and found cached knn results, Pegasus will use cached results and do not recompute.
        * Otherwise, compute kNN irrespective of caching status.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm[rep + "_knn_indices"]``: kNN index matrix. Row i is the index list of kNN of cell i (excluding itself), sorted from nearest to farthest.
        * ``data.obsm[rep + "_knn_distances"]``: kNN distance matrix. Row i is the distance list of kNN of cell i (excluding itselt), sorted from smallest to largest.

    Update ``data.obsp``:
        * ``data.obsp["W_" + rep]``: kNN graph of the data in terms of affinity matrix.

    Examples
    --------
    >>> pg.neighbors(data)
    """

    # calculate kNN
    rep = update_rep(rep)
    indices, distances = get_neighbors(
        data,
        K=K,
        rep=rep,
        n_jobs=n_jobs,
        random_state=random_state,
        full_speed=full_speed,
        use_cache=use_cache,
    )

    # calculate affinity matrix
    W = calculate_affinity_matrix(indices[:, 0 : K - 1], distances[:, 0 : K - 1])
    data.obsp["W_" + rep] = W
    # pop out jump method values
    data.uns.pop(f"{rep}_jump_values", None)
    data.uns.pop(f"{rep}_optimal_k", None)


def calc_kBET_for_one_chunk(knn_indices, attr_values, ideal_dist, K):
    dof = ideal_dist.size - 1

    ns = knn_indices.shape[0]
    results = np.zeros((ns, 2))
    for i in range(ns):
        observed_counts = (
            pd.Series(attr_values[knn_indices[i, :]]).value_counts(sort=False).values
        )
        expected_counts = ideal_dist * K
        stat = np.sum(
            np.divide(
                np.square(np.subtract(observed_counts, expected_counts)),
                expected_counts,
            )
        )
        p_value = 1 - chi2.cdf(stat, dof)
        results[i, 0] = stat
        results[i, 1] = p_value

    return results

def calc_kBET(
    data: MultimodalData,
    attr: str,
    rep: str = "pca",
    K: int = 25,
    alpha: float = 0.05,
    n_jobs: int = -1,
    random_state: int = 0,
    temp_folder: str = None,
    use_cache: bool = True,
) -> Tuple[float, float, float]:
    """Calculate the kBET metric of the data regarding a specific sample attribute and embedding.

    The kBET metric is defined in [Büttner18]_, which measures if cells from different samples mix well in their local neighborhood.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    attr: ``str``
        The sample attribute to consider. Must exist in ``data.obs``.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation to be used. The key ``'X_' + rep`` must exist in ``data.obsm``. By default, use PCA coordinates.

    K: ``int``, optional, default: ``25``
        Number of nearest neighbors, using L2 metric.

    alpha: ``float``, optional, default: ``0.05``
        Acceptance rate threshold. A cell is accepted if its kBET p-value is greater than or equal to ``alpha``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads used. If ``-1``, use all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    temp_folder: ``str``, optional, default: ``None``
        Temporary folder for joblib execution.

    use_cache: ``bool``, optional, default: ``True``
        If use cache results for kNN.

    Returns
    -------
    stat_mean: ``float``
        Mean kBET chi-square statistic over all cells.

    pvalue_mean: ``float``
        Mean kBET p-value over all cells.

    accept_rate: ``float``
        kBET Acceptance rate of the sample.

    Examples
    --------
    >>> pg.calc_kBET(data, attr = 'Channel')

    >>> pg.calc_kBET(data, attr = 'Channel', rep = 'umap')
    """
    assert attr in data.obs
    if data.obs[attr].dtype.name != "category":
        data.obs[attr] = pd.Categorical(data.obs[attr])

    ideal_dist = (
        data.obs[attr].value_counts(normalize=True, sort=False).values
    )  # ideal no batch effect distribution
    nsample = data.shape[0]
    nbatch = ideal_dist.size

    attr_values = data.obs[attr].values.copy()
    attr_values.categories = range(nbatch)

    indices, distances = get_neighbors(
        data, K=K, rep=rep, n_jobs=n_jobs, random_state=random_state, use_cache=use_cache,
    )
    knn_indices = np.concatenate(
        (np.arange(nsample).reshape(-1, 1), indices[:, 0 : K - 1]), axis=1
    )  # add query as 1-nn

    # partition into chunks
    n_jobs = min(eff_n_jobs(n_jobs), nsample)
    starts = np.zeros(n_jobs + 1, dtype=int)
    quotient = nsample // n_jobs
    remainder = nsample % n_jobs
    for i in range(n_jobs):
        starts[i + 1] = starts[i] + quotient + (1 if i < remainder else 0)

    from joblib import Parallel, delayed, parallel_backend
    with parallel_backend("loky", inner_max_num_threads=1):
        kBET_arr = np.concatenate(
            Parallel(n_jobs=n_jobs, temp_folder=temp_folder)(
                delayed(calc_kBET_for_one_chunk)(
                    knn_indices[starts[i] : starts[i + 1], :], attr_values, ideal_dist, K
                )
                for i in range(n_jobs)
            )
        )

    res = kBET_arr.mean(axis=0)
    stat_mean = res[0]
    pvalue_mean = res[1]
    accept_rate = (kBET_arr[:, 1] >= alpha).sum() / nsample

    return (stat_mean, pvalue_mean, accept_rate)


def calc_kSIM(
    data: MultimodalData,
    attr: str,
    rep: str = "pca",
    K: int = 25,
    min_rate: float = 0.9,
    n_jobs: int = -1,
    random_state: int = 0,
    use_cache: bool = True,
) -> Tuple[float, float]:
    """Calculate the kSIM metric of the data regarding a specific sample attribute and embedding.

    The kSIM metric is defined in [Li20]_, which measures if a sample attribute is not diffused too much in each cell's local neighborhood.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    attr: ``str``
        The sample attribute to consider. Must exist in ``data.obs``.

    rep: ``str``, optional, default: ``"pca"``
        The embedding representation to consider. The key ``'X_' + rep`` must exist in ``data.obsm``.

    K: ``int``, optional, default: ``25``
        The number of nearest neighbors to be considered.

    min_rate: ``float``, optional, default: ``0.9``
        Acceptance rate threshold. A cell is accepted if its kSIM rate is larger than or equal to ``min_rate``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads used. If ``-1``, use all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    use_cache: ``bool``, optional, default: ``True``
        If use cache results for kNN.

    Returns
    -------
    kSIM_mean: ``float``
        Mean kSIM rate over all the cells.

    kSIM_accept_rate: ``float``
        kSIM Acceptance rate of the sample.

    Examples
    --------
    >>> pg.calc_kSIM(data, attr = 'cell_type')

    >>> pg.calc_kSIM(data, attr = 'cell_type', rep = 'umap')
    """
    assert attr in data.obs
    nsample = data.shape[0]

    indices, distances = get_neighbors(
        data, K=K, rep=rep, n_jobs=n_jobs, random_state=random_state, use_cache=use_cache,
    )
    knn_indices = np.concatenate(
        (np.arange(nsample).reshape(-1, 1), indices[:, 0 : K - 1]), axis=1
    )  # add query as 1-nn

    labels = np.reshape(data.obs[attr].values[knn_indices.ravel()], (-1, K))
    same_labs = labels == labels[:, 0].reshape(-1, 1)
    correct_rates = same_labs.sum(axis=1) / K

    kSIM_mean = correct_rates.mean()
    kSIM_accept_rate = (correct_rates >= min_rate).sum() / nsample

    return (kSIM_mean, kSIM_accept_rate)
