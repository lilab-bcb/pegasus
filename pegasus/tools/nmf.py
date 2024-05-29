import numpy as np
import pandas as pd

import numba
from numba import njit
from numba.typed import List as numbaList

from typing import List, Union, Tuple
from pegasusio import UnimodalData, MultimodalData
from pegasus.tools import slicing, eff_n_jobs, calculate_nearest_neighbors, check_batch_key

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



def _select_and_scale_features(
    data: Union[MultimodalData, UnimodalData],
    features: str = "highly_variable_features",
    space: str = "log",
    batch: str = None,
) -> Union[np.ndarray, List[np.ndarray]]:
    """ Subset the features as a dense matrix and apply L2 normalization via each column.
    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    features: ``str``, optional, default: ``None``
        a keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.
    space: ``str``, optional, default: ``log``
        Choose from ``log`` or ``expression``. If ``expression``, transfer back to expression space (invert log)
    batch: ``str``, optional, default: ``None``
        a keyword in ``data.obs``, which refers to batches. If batch is not None, return scaled matrix as a list of matrices (one per batch)
    Returns
    -------
    X: ``np.ndarray`` or ``List[np.ndarray]``
        The scaled data matrix or matrices
    """
    if features is not None:
        assert features in data.var
        X = data.X[:, data.var[features].values]
    else:
        X = data.X

    X = slicing(X, copy=True)

    if space != "log":
        np.expm1(X, out = X) # convert back to expression space

    if batch == None:
        X /= np.linalg.norm(X, axis=0)
        return X

    cumsums = pd.concat([pd.Series({'__start': 0}), data.obs[batch].values.value_counts()]).cumsum().values
    Xs = []
    for i in range(cumsums.size - 1):
        X_i = X[cumsums[i]:cumsums[i+1]]
        scale = np.linalg.norm(X_i, axis=0)
        scale[scale == 0.0] = 1.0
        X_i /= scale
        Xs.append(X_i)

    return Xs


@timer(logger=logger)
def nmf(
    data: Union[MultimodalData, UnimodalData],
    n_components: int = 20,
    features: str = "highly_variable_features",
    space: str = "log",
    init: str = "nndsvdar",
    algo: str = "halsvar",
    mode: str = "batch",
    tol: float = 1e-4,
    use_gpu: bool = False,
    alpha_W: float = 0.0,
    l1_ratio_W: float = 0.0,
    alpha_H: float = 0.0,
    l1_ratio_H: float = 0.0,
    fp_precision: str = "float",
    online_chunk_size: int = 5000,
    n_jobs: int = -1,
    random_state: int = 0,
) -> None:
    """Perform Nonnegative Matrix Factorization (NMF) to the data using Frobenius norm. Steps include select features and L2 normalization and NMF and L2 normalization of resulting coordinates.

    The calculation uses `nmf-torch <https://github.com/lilab-bcb/nmf-torch>`_ package.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``50``.
        Number of Principal Components to get.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for nmf.

    max_value: ``float``, optional, default: ``None``.
        The threshold to truncate data symmetrically after scaling. If ``None``, do not truncate.

    space: ``str``, optional, default: ``log``.
        Choose from ``log`` and ``expression``. ``log`` works on log-transformed expression space; ``expression`` works on the original expression space (normalized by total UMIs).

    init: ``str``, optional, default: ``nndsvdar``.
        Method to initialize NMF. Options are 'random', 'nndsvd', 'nndsvda' and 'nndsvdar'.

    algo: ``str``, optional, default: ``halsvar``
        Choose from ``mu`` (Multiplicative Update), ``hals`` (Hierarchical Alternative Least Square), ``halsvar`` (HALS variant, use HALS to mimic ``bpp`` and can get better convergence for sometimes) and ``bpp`` (alternative non-negative least squares with Block Principal Pivoting method).

    mode: ``str``, optional, default: ``batch``
        Learning mode. Choose from ``batch`` and ``online``. Notice that ``online`` only works when ``beta=2.0``. For other beta loss, it switches back to ``batch`` method.

    tol: ``float``, optional, default: ``1e-4``
        The toleration used for convergence check.

    use_gpu: ``bool``, optional, default: ``False``
        If ``True``, use GPU if available. Otherwise, use CPU only.

    alpha_W: ``float``, optional, default: ``0.0``
        A numeric scale factor which multiplies the regularization terms related to W.
        If zero or negative, no regularization regarding W is considered.

    l1_ratio_W: ``float``, optional, default: ``0.0``
        The ratio of L1 penalty on W, must be between 0 and 1. And thus the ratio of L2 penalty on W is (1 - l1_ratio_W).

    alpha_H: ``float``, optional, default: ``0.0``
        A numeric scale factor which multiplies the regularization terms related to H.
        If zero or negative, no regularization regarding H is considered.

    l1_ratio_H: ``float``, optional, default: ``0.0``
        The ratio of L1 penalty on W, must be between 0 and 1. And thus the ratio of L2 penalty on H is (1 - l1_ratio_H).

    fp_precision: ``str``, optional, default: ``float``
        The numeric precision on the results. Choose from ``float`` and ``double``.

    online_chunk_size: ``int``, optional, default: ``int``
        The chunk / mini-batch size for online learning. Only works when ``mode='online'``.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``.
        Random seed to be set for reproducing result.

    Returns
    -------
    ``None``.

    Update ``data.obsm``:

        * ``data.obsm["X_nmf"]``: Scaled NMF coordinates of shape ``(n_cells, n_components)``. Each column has a unit variance.

        * ``data.obsm["H"]``: The coordinate factor matrix of shape ``(n_cells, n_components)``.

    Update ``data.uns``:

        * ``data.uns["W"]``: The feature factor matrix of shape ``(n_HVFs, n_components)``.

        * ``data.uns["nmf_err"]``: The NMF loss.

        * ``data.uns["nmf_features"]``: Record the features used to perform NMF analysis.

    Examples
    --------
    >>> pg.nmf(data)
    """
    X = _select_and_scale_features(data, features=features, space=space)

    try:
        from nmf import run_nmf
    except ImportError as e:
        import sys
        logger.error(f"{e}\nNeed NMF-Torch! Try 'pip install nmf-torch'.")
        sys.exit(-1)

    H, W, err = run_nmf(
        X,
        n_components=n_components,
        init=init,
        algo=algo,
        mode=mode,
        tol=tol,
        n_jobs=eff_n_jobs(n_jobs),
        random_state=random_state,
        use_gpu=use_gpu,
        alpha_W=alpha_W,
        l1_ratio_W=l1_ratio_W,
        alpha_H=alpha_H,
        l1_ratio_H=l1_ratio_H,
        fp_precision=fp_precision,
        online_chunk_size=online_chunk_size,
    )

    data.uns["nmf_features"] = features # record which feature to use
    data.uns["W"] = np.ascontiguousarray(W.T, dtype=np.float32) # cannot be varm because numbers of features are not the same
    data.uns["nmf_err"] = err

    data.obsm["H"] = np.ascontiguousarray(H, dtype=np.float32)
    H = data.obsm["H"]
    data.obsm["X_nmf"] = H / np.linalg.norm(H, axis=0)


@timer(logger=logger)
def find_nmf_programs(
    data: Union[MultimodalData, UnimodalData],
    n_range: Tuple[int, int] = (4, 9),
    n_rep: int = 10,
    features: str = "highly_variable_features",
    space: str = "log",
    init: str = "random",
    algo: str = "halsvar",
    mode: str = "batch",
    tol: float = 1e-4,
    use_gpu: bool = False,
    alpha_W: float = 0.0,
    l1_ratio_W: float = 0.0,
    alpha_H: float = 0.01,
    l1_ratio_H: float = 1.0,
    fp_precision: str = "float",
    online_chunk_size: int = 5000,
    n_jobs: int = -1,
    random_state: int = 0,
) -> Tuple[list, list, list, list]:
    """Perform Nonnegative Matrix Factorization (NMF) to the data using Frobenius norm. Steps include select features and L2 normalization and NMF and L2 normalization of resulting coordinates.

    The calculation uses `nmf-torch <https://github.com/lilab-bcb/nmf-torch>`_ package.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    n_range: ``Tuple[int, int]``, optional, default: ``(4, 9)``.
        Number of ranks to iterate over.

    n_rep: ``int``, optional, default: 10
        Number of reruns for each value in n_range.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for nmf.

    max_value: ``float``, optional, default: ``None``.
        The threshold to truncate data symmetrically after scaling. If ``None``, do not truncate.

    space: ``str``, optional, default: ``log``.
        Choose from ``log`` and ``expression``. ``log`` works on log-transformed expression space; ``expression`` works on the original expression space (normalized by total UMIs).

    init: ``str``, optional, default: ``random``.
        Method to initialize NMF. Options are 'random', 'nndsvd', 'nndsvda' and 'nndsvdar'.

    algo: ``str``, optional, default: ``halsvar``
        Choose from ``mu`` (Multiplicative Update), ``hals`` (Hierarchical Alternative Least Square), ``halsvar`` (HALS variant, use HALS to mimic ``bpp`` and can get better convergence for sometimes) and ``bpp`` (alternative non-negative least squares with Block Principal Pivoting method).

    mode: ``str``, optional, default: ``batch``
        Learning mode. Choose from ``batch`` and ``online``. Notice that ``online`` only works when ``beta=2.0``. For other beta loss, it switches back to ``batch`` method.

    tol: ``float``, optional, default: ``1e-4``
        The toleration used for convergence check.

    use_gpu: ``bool``, optional, default: ``False``
        If ``True``, use GPU if available. Otherwise, use CPU only.

    alpha_W: ``float``, optional, default: ``0.0``
        A numeric scale factor which multiplies the regularization terms related to W.
        If zero or negative, no regularization regarding W is considered.

    l1_ratio_W: ``float``, optional, default: ``0.0``
        The ratio of L1 penalty on W, must be between 0 and 1. And thus the ratio of L2 penalty on W is (1 - l1_ratio_W).

    alpha_H: ``float``, optional, default: ``0.01``
        A numeric scale factor which multiplies the regularization terms related to H.
        If zero or negative, no regularization regarding H is considered.

    l1_ratio_H: ``float``, optional, default: ``1.0``
        The ratio of L1 penalty on W, must be between 0 and 1. And thus the ratio of L2 penalty on H is (1 - l1_ratio_H).

    fp_precision: ``str``, optional, default: ``float``
        The numeric precision on the results. Choose from ``float`` and ``double``.

    online_chunk_size: ``int``, optional, default: ``int``
        The chunk / mini-batch size for online learning. Only works when ``mode='online'``.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``.
        Random seed to be set for reproducing result.

    Returns
    -------
    Hs: best H for each k in n_range
    Ws: best W for each k in n_range
    errs: best err for each k in n_range
    coph_corrs: cophenetic correlation coefficients for each k in n_range

    Examples
    --------
    >>> Hs, Ws, errs, coph_corrs = pg.find_nmf_programs(data)
    """
    X = _select_and_scale_features(data, features=features, space=space)

    try:
        from nmf import run_nmf
        from scipy.cluster.hierarchy import linkage, cophenet
        from scipy.spatial.distance import squareform
    except ImportError as e:
        import sys
        logger.error(f"{e}\nNeed NMF-Torch! Try 'pip install nmf-torch'.")
        sys.exit(-1)

    Hs = []
    Ws = []
    errs = []
    coph_corrs = []

    rng = np.random.default_rng(random_state)
    BIG_NUM = 1000000000
    mats_conn = np.zeros((n_rep, X.shape[0], X.shape[0])) # connectivity matrices

    for k in range(n_range[0], n_range[1] + 1):
        print(f"Begin k={k}:")

        H_best = W_best = None
        err_best = 1e100

        for i in range(n_rep):
            H, W, err = run_nmf(
                X,
                n_components=k,
                init=init,
                algo=algo,
                mode=mode,
                tol=tol,
                n_jobs=eff_n_jobs(n_jobs),
                random_state=rng.integers(BIG_NUM),
                use_gpu=use_gpu,
                alpha_W=alpha_W,
                l1_ratio_W=l1_ratio_W,
                alpha_H=alpha_H,
                l1_ratio_H=l1_ratio_H,
                fp_precision=fp_precision,
                online_chunk_size=online_chunk_size,
            )
            
            if err_best > err:
                err_best = err
                H_best = H
                W_best = W

            clusters = H.argmax(axis=1)
            mats_conn[i] = clusters.reshape((-1, 1)) == clusters.reshape((1, -1))

        consensus = mats_conn.mean(axis=0)
        Y = squareform(1.0 - consensus)
        Z = linkage(Y, method='average')
        coph_corr = cophenet(Z, Y)[0]

        Hs.append(H_best)
        Ws.append(W_best)
        errs.append(err_best)
        coph_corrs.append(coph_corr)

    return Hs, Ws, errs, coph_corrs


@njit(fastmath=True, cache=True)
def _refine_cluster(clusters, indices, ncluster):
    results = np.zeros_like(clusters)
    counter = np.zeros(ncluster, dtype=numba.i4)
    bins = np.zeros(ncluster+1, dtype=numba.i4)
    for i in range(clusters.size):
        counter[:] = 0
        counter[clusters[i]] += 1
        for j in indices[i]:
            counter[clusters[j]] += 1
        clust_id = np.argmax(counter)
        results[i] = clust_id
        bins[clust_id+1] += 1
    return results, np.cumsum(bins)


@njit(fastmath=True, cache=True)
def _quantile_norm(Hs, csums, ids_by_clusts, nbatch, ref_batch, ncluster, min_cells=20, quantiles=50):
    qs = np.linspace(0, 1, quantiles+1) # Generate quantiles

    # Prepare reference batch
    ref_quantiles = []
    ref_sizes = np.zeros(ncluster, dtype=numba.i4)

    Href = Hs[ref_batch]
    csum_ref = csums[ref_batch]
    ids_ref = ids_by_clusts[ref_batch]

    for j in range(ncluster):
        start = csum_ref[j]
        end = csum_ref[j+1]
        ref_sizes[j] = end - start
        if ref_sizes[j] < min_cells:
            ref_quantiles.append([np.zeros(0)])
        else:
            ref_q = []
            Hsub = Href[ids_ref[start:end]]
            for k in range(ncluster):
                values = Hsub[:, k]
                if np.unique(values).size == 1:
                    ref_q.append(np.zeros(0))
                else:
                    ref_q.append(np.quantile(values, qs))
            ref_quantiles.append(ref_q)

    # Quantile normalization
    for i in range(nbatch):
        if i != ref_batch:
            H = Hs[i]
            csum = csums[i]
            ids_by_clust = ids_by_clusts[i]

            for j in range(ncluster):
                start = csum[j]
                end = csum[j+1]
                size = end - start
                if size >= min_cells and ref_sizes[j] >= min_cells:
                    ids = ids_by_clust[start:end]
                    Hsub = H[ids]
                    for k in range(ncluster):
                        if ref_quantiles[j][k].size == 0:
                            H[ids, k] = 0.0
                        else:
                            values = Hsub[:, k]
                            if np.unique(values).size == 1:
                                H[ids, k] = 0.0
                            else:
                                quants = np.quantile(values, qs)
                                H[ids, k] = np.interp(values, quants, ref_quantiles[j][k])


@timer(logger=logger)
def integrative_nmf(
    data: Union[MultimodalData, UnimodalData],
    batch: str = "Channel",
    n_components: int = 20,
    features: str = "highly_variable_features",
    space: str = "log",
    algo: str = "halsvar",
    mode: str = "online",
    tol: float = 1e-4,
    use_gpu: bool = False,
    lam: float = 5.0,
    fp_precision: str = "float",
    online_chunk_size: int = 5000,
    n_jobs: int = -1,
    random_state: int = 0,
    quantile_norm: bool = True,
) -> str:
    """Perform Integrative Nonnegative Matrix Factorization (iNMF) [Yang16]_ for data integration.

    The calculation uses `nmf-torch <https://github.com/lilab-bcb/nmf-torch>`_ .

    This function assumes that cells in each batch are adjacent to each other.
    In addition, it will scale each batch with L2 norm separately. The resulting Hs will also be scaled with L2 norm.
    If ``quantile_norm=True``, quantile normalization will be additionally performed.

    See [Welch19]_ and [Gao21]_ for preprocessing and normalization details.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    batch: ``str``, optional, default: ``"Channel"``.
        Which attribute in data.obs field represents batches, default is "Channel".

    n_components: ``int``, optional, default: ``50``.
        Number of Principal Components to get.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for integrative_nmf.

    space: ``str``, optional, default: ``log``.
        Choose from ``log`` and ``expression``. ``log`` works on log-transformed expression space; ``expression`` works on the original expression space (normalized by total UMIs).

    algo: ``str``, optional, default: ``halsvar``
        Choose from ``mu`` (Multiplicative Update), ``halsvar`` (HALS variant that mimic bpp but faster) and ``bpp`` (alternative non-negative least squares with Block Principal Pivoting method).

    mode: ``str``, optional, default: ``online``
        Learning mode. Choose from ``batch`` and ``online``. Notice that ``online`` only works when ``beta=2.0``. For other beta loss, it switches back to ``batch`` method.

    tol: ``float``, optional, default: ``1e-4``
        The toleration used for convergence check.

    use_gpu: ``bool``, optional, default: ``False``
        If ``True``, use GPU if available. Otherwise, use CPU only.

    lam: ``float``, optional, default: ``5.0``
        The coefficient for regularization terms. If ``0``, then no regularization will be performed.

    fp_precision: ``str``, optional, default: ``float``
        The numeric precision on the results. Choose from ``float`` and ``double``.

    online_chunk_size: ``int``, optional, default: ``5000``
        The chunk / mini-batch size for online learning. Only works when ``mode='online'``.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``.
        Random seed to be set for reproducing result.

    quantile_norm: ``bool``, optioanl, default: ``True``.
        Perform quantile normalization as described in Gao et al. Nature Biotech 2021. Cluster refinement K=20; min_cells=20; quantiles = 50.

    Returns
    -------
    out_rep: ``str``
        The keyword in ``data.obsm`` referring to the embedding calculated by integrative NMF algorithm. out_rep is always equal to "inmf"


    Update ``data.obsm``:

        * ``data.obsm["X_inmf"]``: Scaled and possibly quantile normalized iNMF coordinates.

        * ``data.obsm["H"]``: The concatenation of coordinate factor matrices of shape ``(n_cells, n_components)``.

    Update ``data.uns``:

        * ``data.uns["W"]``: The feature factor matrix of shape ``(n_HVFs, n_components)``.

        * ``data.uns["V"]``: The batch specific feature factor matrices as one tensor of shape ``(n_batches, n_components, n_HVFs)``.

        * ``data.uns["inmf_err"]``: The iNMF loss.

        * ``data.uns["inmf_features"]``: Record the features used to perform iNMF analysis.

    Examples
    --------
    >>> pg.integrative_nmf(data)
    """
    if not check_batch_key(data, batch, "Cannot apply integrative_nmf!"):
        return "pca"

    Xs = _select_and_scale_features(data, features=features, space=space, batch=batch)

    try:
        from nmf import integrative_nmf
    except ImportError as e:
        import sys
        logger.error(f"{e}\nNeed NMF-Torch! Try 'pip install nmf-torch'.")
        sys.exit(-1)

    n_jobs = eff_n_jobs(n_jobs)

    Hs, W, Vs, err = integrative_nmf(
        Xs,
        n_components=n_components,
        algo=algo,
        mode=mode,
        tol=tol,
        n_jobs=n_jobs,
        random_state=random_state,
        use_gpu=use_gpu,
        lam=lam,
        fp_precision=fp_precision,
        online_chunk_size=online_chunk_size,
    )

    # Implementation of algo 3, quantile normalization
    Hs_new = numbaList()
    csums = numbaList()
    ids_by_clusts = numbaList()

    nbatch = len(Hs)
    rg = np.random.default_rng(random_state)
    seeds = rg.integers(4294967295, size=nbatch)
    ref_batch = max_size = -1
    for i in range(nbatch):
        h_norm = np.linalg.norm(Hs[i], axis=0)
        idx_h_zeros = np.where(h_norm==0)[0]
        if idx_h_zeros.size > 0:
            # Set norm 0 to 1 to avoid divide by zero issue
            h_norm[idx_h_zeros] = 1.0
        H_new = np.ascontiguousarray(Hs[i] / h_norm, dtype=np.float32) # Scale H
        Hs_new.append(H_new) # Append scaled H

        if not quantile_norm:
            continue

        clusters = np.argmax(H_new, axis=1) # Assign cluster
        indices, _, _ = calculate_nearest_neighbors(H_new, K=20, n_jobs=n_jobs, random_state=seeds[i]) # KNN with K=20
        clusters, csum = _refine_cluster(clusters, indices, n_components) # Refine cluster
        csums.append(csum)
        ids_by_clusts.append(np.argsort(clusters, kind='stable'))

        if H_new.shape[0] > max_size: # Find ref batch
            max_size = H_new.shape[0]
            ref_batch = i

    if quantile_norm:
        _quantile_norm(Hs_new, csums, ids_by_clusts, nbatch, ref_batch, n_components) # quantile normalization

    data.uns["inmf_features"] = features # record which feature to use
    data.uns["W"] = np.ascontiguousarray(W.T, dtype=np.float32)  # cannot be varm because numbers of features are not the same
    data.uns["V"] = np.array(Vs)
    data.uns["inmf_err"] = err

    data.obsm["H"] = np.concatenate(Hs)
    data.obsm["X_inmf"] = np.concatenate(Hs_new)

    return "inmf"
