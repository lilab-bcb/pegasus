import time
import numpy as np
import scipy
from joblib import effective_n_jobs
from MulticoreTSNE import MulticoreTSNE as TSNE
import umap as umap_module
import forceatlas2 as fa2
import logging

logger = logging.getLogger('sccloud')

from scCloud.tools import (
    neighbors,
    net_train_and_predict,
    calculate_nearest_neighbors,
    calculate_affinity_matrix,
    construct_graph,
)


def calc_tsne(
    X,
    n_jobs,
    n_components,
    perplexity,
    early_exaggeration,
    learning_rate,
    random_state,
    init="random",
    n_iter=1000,
    n_iter_early_exag=250,
):
    """
    TODO: Typing
    """
    tsne = TSNE(
        n_jobs=n_jobs,
        n_components=n_components,
        perplexity=perplexity,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        random_state=random_state,
        verbose=1,
        init=init,
        n_iter=n_iter,
        n_iter_early_exag=n_iter_early_exag,
    )
    X_tsne = tsne.fit_transform(X)
    logger.info("Final error = {}".format(tsne.kl_divergence_))
    return X_tsne


def calc_fitsne(
    X,
    nthreads,
    no_dims,
    perplexity,
    early_exag_coeff,
    learning_rate,
    rand_seed,
    initialization=None,
    max_iter=1000,
    stop_early_exag_iter=250,
    mom_switch_iter=250,
):
    """
    TODO: Typing
    """
    # FItSNE will change X content

    # Check if fftw3 is installed.
    import ctypes.util

    fftw3_loc = ctypes.util.find_library("fftw3")
    if fftw3_loc is None:
        raise Exception("Please install 'fftw3' first to use the FIt-SNE feature!")

    from fitsne import FItSNE

    return FItSNE(
        X.astype("float64"),
        nthreads=nthreads,
        no_dims=no_dims,
        perplexity=perplexity,
        early_exag_coeff=early_exag_coeff,
        learning_rate=learning_rate,
        rand_seed=rand_seed,
        initialization=initialization,
        max_iter=max_iter,
        stop_early_exag_iter=stop_early_exag_iter,
        mom_switch_iter=mom_switch_iter,
    )


# Running umap using our own kNN indices
def calc_umap(
    X,
    n_components,
    n_neighbors,
    min_dist,
    spread,
    random_state,
    init="spectral",
    n_epochs=None,
    learning_rate=1.0,
    knn_indices=None,
    knn_dists=None,
):
    """
    TODO: Typing
    """
    umap_obj = umap_module.UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        random_state=random_state,
        init=init,
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        verbose=True,
    )

    embedding = None
    if X.shape[0] < 4096 or knn_indices is None:
        embedding = umap_obj.fit_transform(X)
        logger.info("using umap kNN graph {}".format(X.shape[0]))
    else:
        assert knn_dists is not None
        # preprocessing codes adopted from UMAP's umap_.py fit function in order to use our own kNN graphs
        from sklearn.utils import check_random_state, check_array

        X = check_array(X, dtype=np.float32, accept_sparse="csr")
        umap_obj._raw_data = X
        if umap_obj.a is None or umap_obj.b is None:
            umap_obj._a, umap_obj._b = umap_module.umap_.find_ab_params(
                umap_obj.spread, umap_obj.min_dist
            )
        else:
            umap_obj._a = umap_obj.a
            umap_obj._b = umap_obj.b
        umap_obj._metric_kwds = (
            umap_obj.metric_kwds if umap_obj.metric_kwds is not None else {}
        )
        umap_obj._target_metric_kwds = {}
        _init = (
            check_array(umap_obj.init, dtype=np.float32, accept_sparse=False)
            if isinstance(umap_obj.init, np.ndarray)
            else umap_obj.init
        )
        umap_obj._initial_alpha = umap_obj.learning_rate
        umap_obj._validate_parameters()

        if umap_obj.verbose:
            logger.info(str(umap_obj))

        if scipy.sparse.isspmatrix_csr(X):
            if not X.has_sorted_indices:
                X.sort_indices()
            umap_obj._sparse_data = True
        else:
            umap_obj._sparse_data = False

        _random_state = check_random_state(umap_obj.random_state)

        if umap_obj.verbose:
            logger.info("Construct fuzzy simplicial set")

        umap_obj._small_data = False
        umap_obj.graph_ = umap_module.umap_.fuzzy_simplicial_set(
            X,
            umap_obj.n_neighbors,
            _random_state,
            umap_obj.metric,
            umap_obj._metric_kwds,
            knn_indices,
            knn_dists,
            umap_obj.angular_rp_forest,
            umap_obj.set_op_mix_ratio,
            umap_obj.local_connectivity,
            umap_obj.verbose,
        )

        _n_epochs = umap_obj.n_epochs if umap_obj.n_epochs is not None else 0
        if umap_obj.verbose:
            logger.info("Construct embedding")
        embedding = umap_module.umap_.simplicial_set_embedding(
            X,
            umap_obj.graph_,
            umap_obj.n_components,
            umap_obj._initial_alpha,
            umap_obj._a,
            umap_obj._b,
            umap_obj.repulsion_strength,
            umap_obj.negative_sample_rate,
            _n_epochs,
            _init,
            _random_state,
            umap_obj.metric,
            umap_obj._metric_kwds,
            umap_obj.verbose,
        )

    return embedding


def calc_force_directed_layout(
    W,
    file_name,
    n_jobs,
    target_change_per_node,
    target_steps,
    is3d,
    memory,
    random_state,
    init=None,
):
    """
    TODO: Typing
    """
    G = construct_graph(W)
    return fa2.forceatlas2(file_name, graph=G, n_jobs=n_jobs, target_change_per_node=target_change_per_node, \
        target_steps=target_steps, is3d=is3d, memory=memory, random_state=random_state, init=init)


def tsne(
    data: "AnnData",
    rep: "str" = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    perplexity: float = 30,
    early_exaggeration: int = 12,
    learning_rate: float = 1000,
    random_state: int = 0,
    out_basis: str = "tsne",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()

    rep_key = "X_" + rep
    if rep_key not in data.obsm.keys():
        raise ValueError("Cannot find {0} matrix. Please run {0} first!".format(rep))

    n_jobs = effective_n_jobs(n_jobs)

    data.obsm["X_" + out_basis] = calc_tsne(
        data.obsm[rep_key],
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
    )

    end = time.time()
    logger.info("t-SNE is calculated. Time spent = {:.2f}s.".format(end - start))


def fitsne(
    data: "AnnData",
    rep: "str" = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    perplexity: float = 30,
    early_exaggeration: int = 12,
    learning_rate: float = 1000,
    random_state: int = 0,
    out_basis: str = "fitsne",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()

    rep_key = "X_" + rep if rep != "CITE-Seq" else rep
    if rep_key not in data.obsm.keys():
        raise ValueError("Cannot find {0} matrix. Please run {0} first!".format(rep))

    n_jobs = effective_n_jobs(n_jobs)

    data.obsm["X_" + out_basis] = calc_fitsne(
        data.obsm[rep_key],
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
    )

    end = time.time()
    logger.info("FIt-SNE is calculated. Time spent = {:.2f}s.".format(end - start))


def umap(
    data: "AnnData",
    rep: "str" = "pca",
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.5,
    spread: float = 1.0,
    random_state: int = 0,
    out_basis: str = "umap",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()

    rep_key = "X_" + rep
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if (
        (indices_key not in data.uns)
        or (distances_key not in data.uns)
        or (n_neighbors > data.uns[indices_key].shape[1] + 1)
    ):
        raise ValueError("Please run neighbors first!")
    if rep_key not in data.obsm.keys():
        raise ValueError("Please run {} first!".format(rep))

    knn_indices = np.insert(
        data.uns[indices_key][:, 0: n_neighbors - 1], 0, range(data.shape[0]), axis=1
    )
    knn_dists = np.insert(
        data.uns[distances_key][:, 0: n_neighbors - 1], 0, 0.0, axis=1
    )
    data.obsm["X_" + out_basis] = calc_umap(
        data.obsm[rep_key],
        n_components,
        n_neighbors,
        min_dist,
        spread,
        random_state,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )

    end = time.time()
    logger.info("UMAP is calculated. Time spent = {:.2f}s.".format(end - start))


def fle(
    data: "AnnData",
    file_name: str,
    n_jobs: int = -1,
    K: int = 50,
    full_speed: bool = False,
    target_change_per_node: float = 2.0,
    target_steps: int = 5000,
    is3d: bool = False,
    memory: int = 8,
    random_state: int = 0,
    out_basis: str = "fle",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()
    n_jobs = effective_n_jobs(n_jobs)

    rep = "diffmap"
    rep_key = "X_" + rep
    W_rep = "W_" + rep

    if rep_key not in data.obsm.keys():
        raise ValueError("Please run diffmap first!")

    if W_rep not in data.uns:
        neighbors(
            data,
            K=K,
            rep=rep,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )

    data.obsm["X_" + out_basis] = calc_force_directed_layout(
        data.uns[W_rep],
        file_name,
        n_jobs,
        target_change_per_node,
        target_steps,
        is3d,
        memory,
        random_state,
    )

    end = time.time()
    logger.info(
        "Force-directed layout is calculated. Time spent = {:.2f}s.".format(end - start)
    )


def select_cells(distances, frac, K=25, alpha=1.0, random_state=0):
    """
    TODO: documentation (not user API)
    """

    start_time = time.time()

    nsample = distances.shape[0]

    if K > distances.shape[1]:
        logger.info(
            "Warning: in select_cells, K = {} > the number of calculated nearest neighbors!\nSet K to {}".format(
                K, distances.shape[1]
            )
        )
        K = distances.shape[1]

    probs = np.zeros(nsample)
    if alpha == 0.0:
        probs[:] = 1.0  # uniform
    elif alpha == 1.0:
        probs[:] = distances[:, K - 1]
    else:
        probs[:] = distances[:, K - 1] ** alpha
    probs /= probs.sum()

    np.random.seed(random_state)
    selected = np.zeros(nsample, dtype=bool)
    selected[
        np.random.choice(nsample, size=int(nsample * frac), replace=False, p=probs)
    ] = True

    end_time = time.time()
    logger.info("select_cells finished. Time spent = {:.2}s.".format(end_time - start_time))

    return selected


def net_tsne(
    data: "AnnData",
    rep: "str" = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    perplexity: float = 30,
    early_exaggeration: int = 12,
    learning_rate: float = 1000,
    random_state: int = 0,
    select_frac: float = 0.1,
    select_K: int = 25,
    select_alpha: float = 1.0,
    net_alpha: float = 0.1,
    polish_learning_frac: float = 0.33,
    polish_n_iter: int = 150,
    out_basis: str = "net_tsne",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()

    rep_key = "X_" + rep
    distances_key = rep + "_knn_distances"

    if rep_key not in data.obsm.keys():
        raise ValueError("Please run {} first!".format(rep))

    if distances_key not in data.uns:
        raise ValueError("Please run neighbors first!")

    n_jobs = effective_n_jobs(n_jobs)

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )

    X = data.obsm[rep_key][selected, :]
    X_tsne = calc_tsne(
        X,
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
    )

    data.uns["X_" + out_basis + "_small"] = X_tsne
    data.obs["ds_selected"] = selected

    Y_init = np.zeros((data.shape[0], 2), dtype=np.float64)
    Y_init[selected, :] = X_tsne
    Y_init[~selected, :] = net_train_and_predict(
        X,
        X_tsne,
        data.obsm[rep_key][~selected, :],
        net_alpha,
        random_state,
        verbose=True,
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    polish_learning_rate = polish_learning_frac * data.shape[0]
    data.obsm["X_" + out_basis] = calc_tsne(
        data.obsm[rep_key],
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        polish_learning_rate,
        random_state,
        init=Y_init,
        n_iter=polish_n_iter,
        n_iter_early_exag=0,
    )

    end = time.time()
    logger.info("Net tSNE is calculated. Time spent = {:.2f}s.".format(end - start))


def net_fitsne(
    data: "AnnData",
    rep: "str" = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    perplexity: float = 30,
    early_exaggeration: int = 12,
    learning_rate: float = 1000,
    random_state: int = 0,
    select_frac: float = 0.1,
    select_K: int = 25,
    select_alpha: float = 1.0,
    net_alpha: float = 0.1,
    polish_learning_frac: float = 0.5,
    polish_n_iter: int = 150,
    out_basis: "str" = "net_fitsne",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()

    rep_key = "X_" + rep
    distances_key = rep + "_knn_distances"

    if rep_key not in data.obsm.keys():
        raise ValueError("Please run {} first!".format(rep))

    if distances_key not in data.uns:
        raise ValueError("Please run neighbors first!")

    n_jobs = effective_n_jobs(n_jobs)

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )

    X = data.obsm[rep_key][selected, :]
    X_fitsne = calc_fitsne(
        X,
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
    )

    data.uns["X_" + out_basis + "_small"] = X_fitsne
    data.obs["ds_selected"] = selected

    Y_init = np.zeros((data.shape[0], 2), dtype=np.float64)
    Y_init[selected, :] = X_fitsne
    Y_init[~selected, :] = net_train_and_predict(
        X,
        X_fitsne,
        data.obsm[rep_key][~selected, :],
        net_alpha,
        random_state,
        verbose=True,
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    polish_learning_rate = polish_learning_frac * data.shape[0]
    data.obsm["X_" + out_basis] = calc_fitsne(
        data.obsm[rep_key],
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        polish_learning_rate,
        random_state,
        initialization=Y_init,
        max_iter=polish_n_iter,
        stop_early_exag_iter=0,
        mom_switch_iter=0,
    )

    end = time.time()
    logger.info("Net FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))


def net_umap(
    data: "AnnData",
    rep: "str" = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    spread: float = 1.0,
    random_state: int = 0,
    select_frac: float = 0.1,
    select_K: int = 25,
    select_alpha: float = 1.0,
    full_speed: bool = False,
    net_alpha: float = 0.1,
    polish_learning_rate: float = 10.0,
    polish_n_epochs: int = 30,
    out_basis: str = "net_umap",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()

    rep_key = "X_" + rep
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if rep_key not in data.obsm.keys():
        raise ValueError("Please run {} first!".format(rep))

    if (indices_key not in data.uns) or (distances_key not in data.uns):
        raise ValueError("Please run neighbors first!")

    n_jobs = effective_n_jobs(n_jobs)

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )

    X = data.obsm[rep_key][selected, :]

    ds_indices_key = "ds_" + rep + "_knn_indices"
    ds_distances_key = "ds_" + rep + "_knn_distances"
    if (ds_indices_key not in data.uns) or (ds_distances_key not in data.uns):
        indices, distances = calculate_nearest_neighbors(
            X,
            K=n_neighbors,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )
        data.uns[ds_indices_key] = indices
        data.uns[ds_distances_key] = distances

    knn_indices = np.insert(
        data.uns[ds_indices_key][:, 0: n_neighbors - 1], 0, range(X.shape[0]), axis=1
    )
    knn_dists = np.insert(
        data.uns[ds_distances_key][:, 0: n_neighbors - 1], 0, 0.0, axis=1
    )

    X_umap = calc_umap(
        X,
        n_components,
        n_neighbors,
        min_dist,
        spread,
        random_state,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )

    data.uns["X_" + out_basis + "_small"] = X_umap
    data.obs["ds_selected"] = selected

    Y_init = np.zeros((data.shape[0], 2), dtype=np.float64)
    Y_init[selected, :] = X_umap
    Y_init[~selected, :] = net_train_and_predict(
        X,
        X_umap,
        data.obsm[rep_key][~selected, :],
        net_alpha,
        random_state,
        verbose=True,
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    knn_indices = np.insert(
        data.uns[indices_key][:, 0: n_neighbors - 1], 0, range(data.shape[0]), axis=1
    )
    knn_dists = np.insert(
        data.uns[distances_key][:, 0: n_neighbors - 1], 0, 0.0, axis=1
    )

    data.obsm["X_" + out_basis] = calc_umap(
        data.obsm[rep_key],
        n_components,
        n_neighbors,
        min_dist,
        spread,
        random_state,
        init=Y_init,
        n_epochs=polish_n_epochs,
        learning_rate=polish_learning_rate,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )

    end = time.time()
    logger.info("Net UMAP is calculated. Time spent = {:.2f}s.".format(end - start))


def net_fle(
    data: "AnnData",
    file_name: str,
    n_jobs: int = -1,
    K: int = 50,
    full_speed: bool = False,
    target_change_per_node: float = 2.0,
    target_steps: int = 5000,
    is3d: bool = False,
    memory: int = 8,
    random_state: int = 0,
    select_frac: float = 0.1,
    select_K: int = 25,
    select_alpha: float = 1.0,
    net_alpha: float = 0.1,
    polish_target_steps: int = 1500,
    out_basis: str = "net_fle",
) -> None:
    """
    TODO: Documentation.
    """
    start = time.time()
    n_jobs = effective_n_jobs(n_jobs)

    rep = "diffmap"
    rep_key = "X_" + rep
    W_rep = "W_" + rep

    if rep_key not in data.obsm.keys():
        raise ValueError("Please run diffmap first!")
    if W_rep not in data.uns:
        neighbors(
            data,
            K=K,
            rep=rep,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )

    distances_key = rep + "_knn_distances"
    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )

    X = data.obsm[rep_key][selected, :]

    ds_indices_key = "ds_" + rep + "_knn_indices"
    ds_distances_key = "ds_" + rep + "_knn_distances"
    if (ds_indices_key not in data.uns) or (ds_distances_key not in data.uns):
        indices, distances = calculate_nearest_neighbors(
            X, K=K, n_jobs=n_jobs, random_state=random_state, full_speed=full_speed
        )
        data.uns[ds_indices_key] = indices
        data.uns[ds_distances_key] = distances
    else:
        indices = data.uns[ds_indices_key]
        distances = data.uns[ds_distances_key]

    W = calculate_affinity_matrix(indices, distances)

    X_fle = calc_force_directed_layout(
        W,
        file_name + ".small",
        n_jobs,
        target_change_per_node,
        target_steps,
        is3d,
        memory,
        random_state,
    )

    data.uns["X_" + out_basis + "_small"] = X_fle
    data.obs["ds_diffmap_selected"] = selected

    Y_init = np.zeros((data.shape[0], 2), dtype=np.float64)
    Y_init[selected, :] = X_fle
    Y_init[~selected, :] = net_train_and_predict(
        X,
        X_fle,
        data.obsm[rep_key][~selected, :],
        net_alpha,
        random_state,
        verbose=True,
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    data.obsm["X_" + out_basis] = calc_force_directed_layout(
        data.uns[W_rep],
        file_name,
        n_jobs,
        target_change_per_node,
        polish_target_steps,
        is3d,
        memory,
        random_state,
        init=Y_init,
    )

    end = time.time()
    logger.info("Net FLE is calculated. Time spent = {:.2f}s.".format(end - start))
