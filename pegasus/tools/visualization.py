import time
import numpy as np
import scipy
import umap as umap_module
import forceatlas2 as fa2
import uuid

from pegasusio import MultimodalData
from joblib import effective_n_jobs
try:
    from MulticoreTSNE import MulticoreTSNE as TSNE
except ImportError:
    print("Need Multicore-TSNE!")

from pegasus.tools import (
    update_rep,
    X_from_rep,
    W_from_rep,
    knn_is_cached,
    neighbors,
    net_train_and_predict,
    calculate_nearest_neighbors,
    calculate_affinity_matrix,
    construct_graph,
)

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



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
        umap_obj.graph_, umap_obj._sigmas, umap_obj._rhos = umap_module.umap_.fuzzy_simplicial_set(
            X=X,
            n_neighbors=umap_obj.n_neighbors,
            random_state=_random_state,
            metric=umap_obj.metric,
            metric_kwds=umap_obj._metric_kwds,
            knn_indices=knn_indices,
            knn_dists=knn_dists,
            angular=umap_obj.angular_rp_forest,
            set_op_mix_ratio=umap_obj.set_op_mix_ratio,
            local_connectivity=umap_obj.local_connectivity,
            verbose=umap_obj.verbose,
        )

        _n_epochs = umap_obj.n_epochs if umap_obj.n_epochs is not None else 0
        if umap_obj.verbose:
            logger.info("Construct embedding")
        embedding = umap_module.umap_.simplicial_set_embedding(
            data=X,
            graph=umap_obj.graph_,
            n_components=umap_obj.n_components,
            initial_alpha=umap_obj._initial_alpha,
            a=umap_obj._a,
            b=umap_obj._b,
            gamma=umap_obj.repulsion_strength,
            negative_sample_rate=umap_obj.negative_sample_rate,
            n_epochs=_n_epochs,
            init=_init,
            random_state=_random_state,
            metric=umap_obj.metric,
            metric_kwds=umap_obj._metric_kwds,
            verbose=umap_obj.verbose,
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
    return fa2.forceatlas2(
        file_name,
        graph=G,
        n_jobs=n_jobs,
        target_change_per_node=target_change_per_node,
        target_steps=target_steps,
        is3d=is3d,
        memory=memory,
        random_state=random_state,
        init=init,
    )


@timer(logger=logger)
def tsne(
    data: MultimodalData,
    rep: str = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    perplexity: float = 30,
    early_exaggeration: int = 12,
    learning_rate: float = 1000,
    random_state: int = 0,
    out_basis: str = "tsne",
) -> None:
    """Calculate tSNE embedding using MulticoreTSNE_ package.

    .. _MulticoreTSNE: https://github.com/DmitryUlyanov/Multicore-TSNE

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated tSNE coordinates. By default, generate 2-dimensional data for 2D visualization.

    perplexity: ``float``, optional, default: ``30``
        The perplexity is related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require a larger perplexity.

    early_exaggeration: ``int``, optional, default: ``12``
        Controls how tight natural clusters in the original space are in the embedded space, and how much space will be between them.

    learning_rate: ``float``, optional, default: ``1000``
        The learning rate can be a critical parameter, which should be between 100 and 1000.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    out_basis: ``str``, optional, default: ``"tsne"``
        Key name for calculated tSNE coordinates to store.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: tSNE coordinates of the data.

    Examples
    --------
    >>> pg.tsne(data)
    """

    rep = update_rep(rep)
    n_jobs = effective_n_jobs(n_jobs)

    data.obsm["X_" + out_basis] = calc_tsne(
        X_from_rep(data, rep),
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
    )


@timer(logger=logger)
def fitsne(
    data: MultimodalData,
    rep: str = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    perplexity: float = 30,
    early_exaggeration: int = 12,
    learning_rate: float = 1000,
    random_state: int = 0,
    out_basis: str = "fitsne",
) -> None:
    """Calculate FIt-SNE embedding using fitsne_ package.

    .. _fitsne: https://github.com/KlugerLab/FIt-SNE

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated FI-tSNE coordinates. By default, generate 2-dimensional data for 2D visualization.

    perplexity: ``float``, optional, default: ``30``
        The perplexity is related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require a larger perplexity.

    early_exaggeration: ``int``, optional, default: ``12``
        Controls how tight natural clusters in the original space are in the embedded space, and how much space will be between them.

    learning_rate: ``float``, optional, default: ``1000``
        The learning rate can be a critical parameter, which should be between 100 and 1000.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    out_basis: ``str``, optional, default: ``"fitsne"``
        Key name for calculated FI-tSNE coordinates to store.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: FI-tSNE coordinates of the data.

    Examples
    --------
    >>> pg.fitsne(data)
    """

    rep = update_rep(rep)
    n_jobs = effective_n_jobs(n_jobs)

    data.obsm["X_" + out_basis] = calc_fitsne(
        X_from_rep(data, rep),
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
    )


@timer(logger=logger)
def umap(
    data: MultimodalData,
    rep: str = "pca",
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.5,
    spread: float = 1.0,
    random_state: int = 0,
    out_basis: str = "umap",
) -> None:
    """Calculate UMAP embedding using umap-learn_ package.

    .. _umap-learn: https://github.com/lmcinnes/umap

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated UMAP coordinates. By default, generate 2-dimensional data for 2D visualization.

    n_neighbors: ``int``, optional, default: ``15``
        Number of nearest neighbors considered during the computation.

    min_dist: ``float``, optional, default: ``0.5``
        The effective minimum distance between embedded data points.

    spread: ``float``, optional, default: ``1.0``
        The effective scale of embedded data points.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    out_basis: ``str``, optional, default: ``"umap"``
        Key name for calculated UMAP coordinates to store.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: UMAP coordinates of the data.

    Examples
    --------
    >>> pg.umap(data)
    """
    start = time.time()

    rep = update_rep(rep)
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    X = X_from_rep(data, rep)
    if not knn_is_cached(data, indices_key, distances_key, n_neighbors):
        raise ValueError("Please run neighbors first!")

    knn_indices = np.insert(
        data.uns[indices_key][:, 0 : n_neighbors - 1], 0, range(data.shape[0]), axis=1
    )
    knn_dists = np.insert(
        data.uns[distances_key][:, 0 : n_neighbors - 1], 0, 0.0, axis=1
    )
    data.obsm["X_" + out_basis] = calc_umap(
        X,
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


@timer(logger=logger)
def fle(
    data: MultimodalData,
    file_name: str = None,
    n_jobs: int = -1,
    rep: str = "diffmap",
    K: int = 50,
    full_speed: bool = False,
    target_change_per_node: float = 2.0,
    target_steps: int = 5000,
    is3d: bool = False,
    memory: int = 8,
    random_state: int = 0,
    out_basis: str = "fle",
) -> None:
    """Construct the Force-directed (FLE) graph using ForceAtlas2_ implementation, with Python wrapper as forceatlas2-python_.

    .. _ForceAtlas2: https://github.com/klarman-cell-observatory/forceatlas2
    .. _forceatlas2-python: https://github.com/klarman-cell-observatory/forceatlas2-python

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    file_name: ``str``, optional, default: ``None``
        Temporary file to store the coordinates as the input to forceatlas2. If ``None``, use ``tempfile.mkstemp`` to generate file name.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    rep: ``str``, optional, default: ``"diffmap"``
        Representation of data used for the calculation. By default, use Diffusion Map coordinates. If ``None``, use the count matrix ``data.X``.

    K: ``int``, optional, default: ``50``
        Number of nearest neighbors to be considered during the computation.

    full_speed: ``bool``, optional, default: ``False``
        * If ``True``, use multiple threads in constructing ``hnsw`` index. However, the kNN results are not reproducible.
        * Otherwise, use only one thread to make sure results are reproducible.

    target_change_per_node: ``float``, optional, default: ``2.0``
        Target change per node to stop ForceAtlas2.

    target_steps: ``int``, optional, default: ``5000``
        Maximum number of iterations before stopping the ForceAtlas2 algorithm.

    is3d: ``bool``, optional, default: ``False``
        If ``True``, calculate 3D force-directed layout.

    memory: ``int``, optional, default: ``8``
        Memory size in GB for the Java FA2 component. By default, use 8GB memory.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    out_basis: ``str``, optional, default: ``"fle"``
        Key name for calculated FLE coordinates to store.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: FLE coordinates of the data.

    Examples
    --------
    >>> pg.fle(data)
    """

    if file_name is None:
        import tempfile

        _, file_name = tempfile.mkstemp()

    n_jobs = effective_n_jobs(n_jobs)
    rep = update_rep(rep)

    if ("W_" + rep) not in data.uns:
        neighbors(
            data,
            K=K,
            rep=rep,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )

    data.obsm["X_" + out_basis] = calc_force_directed_layout(
        W_from_rep(data, rep),
        file_name,
        n_jobs,
        target_change_per_node,
        target_steps,
        is3d,
        memory,
        random_state,
    )


@timer(logger=logger)
def select_cells(distances, frac, K=25, alpha=1.0, random_state=0):
    """
    TODO: documentation (not user API)
    """

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

    return selected


@timer(logger=logger)
def net_tsne(
    data: MultimodalData,
    rep: str = "pca",
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
    """Calculate approximated tSNE embedding using Deep Learning model to improve the speed.

    In specific, the deep model used is MLPRegressor_, the *scikit-learn* implementation of Multi-layer Perceptron regressor.

    .. _MLPRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells (``n_obs``) and columns for genes (``n_feature``).

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated tSNE coordinates. By default, generate 2-dimensional data for 2D visualization.

    perplexity: ``float``, optional, default: ``30``
        The perplexity is related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require a larger perplexity.

    early_exaggeration: ``int``, optional, default: ``12``
        Controls how tight natural clusters in the original space are in the embedded space, and how much space will be between them.

    learning_rate: ``float``, optional, default: ``1000``
        The learning rate can be a critical parameter, which should be between 100 and 1000.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    select_frac: ``float``, optional, default: ``0.1``
        Down sampling fraction on the cells.

    select_K: ``int``, optional, default: ``25``
        Number of neighbors to be used to estimate local density for each data point for down sampling.

    select_alpha: ``float``, optional, default: ``1.0``
        Weight the down sample to be proportional to ``radius ** select_alpha``.

    net_alpha: ``float``, optional, default: ``0.1``
        L2 penalty (regularization term) parameter of the deep regressor.

    polish_learning_frac: ``float``, optional, default: ``0.33``
        After running the deep regressor to predict new coordinates, use ``polish_learning_frac`` * ``n_obs`` as the learning rate to polish the coordinates.

    polish_n_iter: ``int``, optional, default: ``150``
        Number of iterations for polishing tSNE run.

    out_basis: ``str``, optional, default: ``"net_tsne"``
        Key name for the approximated tSNE coordinates calculated.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: Net tSNE coordinates of the data.

    Update ``data.obs``:
        * ``data.obs['ds_selected']``: Boolean array to indicate which cells are selected during the down sampling phase.

    Examples
    --------
    >>> pg.net_tsne(data)
    """

    rep = update_rep(rep)
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if not knn_is_cached(data, indices_key, distances_key, select_K):
        raise ValueError("Please run neighbors first!")

    n_jobs = effective_n_jobs(n_jobs)

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )

    X_full = X_from_rep(data, rep)
    X = X_full[selected, :]
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
        X, X_tsne, X_full[~selected, :], net_alpha, random_state, verbose=True
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    polish_learning_rate = polish_learning_frac * data.shape[0]
    data.obsm["X_" + out_basis] = calc_tsne(
        X_full,
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


@timer(logger=logger)
def net_fitsne(
    data: MultimodalData,
    rep: str = "pca",
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
    """Calculate approximated FI-tSNE embedding using Deep Learning model to improve the speed.

    In specific, the deep model used is MLPRegressor_, the *scikit-learn* implementation of Multi-layer Perceptron regressor.

    .. _MLPRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells (``n_obs``) and columns for genes (``n_feature``).

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated tSNE coordinates. By default, generate 2-dimensional data for 2D visualization.

    perplexity: ``float``, optional, default: ``30``
        The perplexity is related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require a larger perplexity.

    early_exaggeration: ``int``, optional, default: ``12``
        Controls how tight natural clusters in the original space are in the embedded space, and how much space will be between them.

    learning_rate: ``float``, optional, default: ``1000``
        The learning rate can be a critical parameter, which should be between 100 and 1000.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    select_frac: ``float``, optional, default: ``0.1``
        Down sampling fraction on the cells.

    select_K: ``int``, optional, default: ``25``
        Number of neighbors to be used to estimate local density for each data point for down sampling.

    select_alpha: ``float``, optional, default: ``1.0``
        Weight the down sample to be proportional to ``radius ** select_alpha``.

    net_alpha: ``float``, optional, default: ``0.1``
        L2 penalty (regularization term) parameter of the deep regressor.

    polish_learning_frac: ``float``, optional, default: ``0.5``
        After running the deep regressor to predict new coordinates, use ``polish_learning_frac`` * ``n_obs`` as the learning rate to polish the coordinates.

    polish_n_iter: ``int``, optional, default: ``150``
        Number of iterations for polishing FI-tSNE run.

    out_basis: ``str``, optional, default: ``"net_fitsne"``
        Key name for the approximated FI-tSNE coordinates calculated.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: Net FI-tSNE coordinates of the data.

    Update ``data.obs``:
        * ``data.obs['ds_selected']``: Boolean array to indicate which cells are selected during the down sampling phase.

    Examples
    --------
    >>> pg.net_fitsne(data)
    """

    rep = update_rep(rep)
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if not knn_is_cached(data, indices_key, distances_key, select_K):
        raise ValueError("Please run neighbors first!")

    n_jobs = effective_n_jobs(n_jobs)

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )
    X_full = X_from_rep(data, rep)
    X = X_full[selected, :]
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
        X, X_fitsne, X_full[~selected, :], net_alpha, random_state, verbose=True
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    polish_learning_rate = polish_learning_frac * data.shape[0]
    data.obsm["X_" + out_basis] = calc_fitsne(
        X_full,
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


@timer(logger=logger)
def net_umap(
    data: MultimodalData,
    rep: str = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.5,
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
    """Calculate approximated UMAP embedding using Deep Learning model to improve the speed.

    In specific, the deep model used is MLPRegressor_, the *scikit-learn* implementation of Multi-layer Perceptron regressor.

    .. _MLPRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated UMAP coordinates. By default, generate 2-dimensional data for 2D visualization.

    n_neighbors: ``int``, optional, default: ``15``
        Number of nearest neighbors considered during the computation.

    min_dist: ``float``, optional, default: ``0.5``
        The effective minimum distance between embedded data points.

    spread: ``float``, optional, default: ``1.0``
        The effective scale of embedded data points.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    select_frac: ``float``, optional, default: ``0.1``
        Down sampling fraction on the cells.

    select_K: ``int``, optional, default: ``25``
        Number of neighbors to be used to estimate local density for each data point for down sampling.

    select_alpha: ``float``, optional, default: ``1.0``
        Weight the down sample to be proportional to ``radius ** select_alpha``.

    full_speed: ``bool``, optional, default: ``False``
        * If ``True``, use multiple threads in constructing ``hnsw`` index. However, the kNN results are not reproducible.
        * Otherwise, use only one thread to make sure results are reproducible.

    net_alpha: ``float``, optional, default: ``0.1``
        L2 penalty (regularization term) parameter of the deep regressor.

    polish_learning_frac: ``float``, optional, default: ``10.0``
        After running the deep regressor to predict new coordinates, use ``polish_learning_frac`` * ``n_obs`` as the learning rate to polish the coordinates.

    polish_n_iter: ``int``, optional, default: ``30``
        Number of iterations for polishing UMAP run.

    out_basis: ``str``, optional, default: ``"net_umap"``
        Key name for calculated UMAP coordinates to store.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: Net UMAP coordinates of the data.

    Update ``data.obs``:
        * ``data.obs['ds_selected']``: Boolean array to indicate which cells are selected during the down sampling phase.

    Examples
    --------
    >>> pg.net_umap(data)
    """

    rep = update_rep(rep)
    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if not knn_is_cached(data, indices_key, distances_key, select_K):
        raise ValueError("Please run neighbors first!")

    n_jobs = effective_n_jobs(n_jobs)

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )
    X_full = X_from_rep(data, rep)
    X = X_full[selected, :]

    ds_indices_key = "ds_" + rep + "_knn_indices"  # ds refers to down-sampling
    ds_distances_key = "ds_" + rep + "_knn_distances"
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
        data.uns[ds_indices_key][:, 0 : n_neighbors - 1], 0, range(X.shape[0]), axis=1
    )
    knn_dists = np.insert(
        data.uns[ds_distances_key][:, 0 : n_neighbors - 1], 0, 0.0, axis=1
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
        X, X_umap, X_full[~selected, :], net_alpha, random_state, verbose=True
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    knn_indices = np.insert(
        data.uns[indices_key][:, 0 : n_neighbors - 1], 0, range(data.shape[0]), axis=1
    )
    knn_dists = np.insert(
        data.uns[distances_key][:, 0 : n_neighbors - 1], 0, 0.0, axis=1
    )

    data.obsm["X_" + out_basis] = calc_umap(
        X_full,
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


@timer(logger=logger)
def net_fle(
    data: MultimodalData,
    file_name: str = None,
    n_jobs: int = -1,
    rep: str = "diffmap",
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
    """Construct the approximated Force-directed (FLE) graph using Deep Learning model to improve the speed.

    In specific, the deep model used is MLPRegressor_, the *scikit-learn* implementation of Multi-layer Perceptron regressor.

    .. _MLPRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    file_name: ``str``, optional, default: ``None``
        Temporary file to store the coordinates as the input to forceatlas2. If ``None``, use ``tempfile.mkstemp`` to generate file name.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    rep: ``str``, optional, default: ``"diffmap"``
        Representation of data used for the calculation. By default, use Diffusion Map coordinates. If ``None``, use the count matrix ``data.X``.

    K: ``int``, optional, default: ``50``
        Number of nearest neighbors to be considered during the computation.

    full_speed: ``bool``, optional, default: ``False``
        * If ``True``, use multiple threads in constructing ``hnsw`` index. However, the kNN results are not reproducible.
        * Otherwise, use only one thread to make sure results are reproducible.

    target_change_per_node: ``float``, optional, default: ``2.0``
        Target change per node to stop ForceAtlas2.

    target_steps: ``int``, optional, default: ``5000``
        Maximum number of iterations before stopping the ForceAtlas2 algorithm.

    is3d: ``bool``, optional, default: ``False``
        If ``True``, calculate 3D force-directed layout.

    memory: ``int``, optional, default: ``8``
        Memory size in GB for the Java FA2 component. By default, use 8GB memory.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    select_frac: ``float``, optional, default: ``0.1``
        Down sampling fraction on the cells.

    select_K: ``int``, optional, default: ``25``
        Number of neighbors to be used to estimate local density for each data point for down sampling.

    select_alpha: ``float``, optional, default: ``1.0``
        Weight the down sample to be proportional to ``radius ** select_alpha``.

    net_alpha: ``float``, optional, default: ``0.1``
        L2 penalty (regularization term) parameter of the deep regressor.

    polish_target_steps: ``int``, optional, default: ``1500``
        After running the deep regressor to predict new coordinate, Number of ForceAtlas2 iterations.

    out_basis: ``str``, optional, default: ``"net_fle"``
        Key name for calculated FLE coordinates to store.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm['X_' + out_basis]``: Net FLE coordinates of the data.

    Update ``data.obs``:
        * ``data.obs['ds_selected']``: Boolean array to indicate which cells are selected during the down sampling phase.

    Examples
    --------
    >>> pg.net_fle(data)
    """

    if file_name is None:
        if file_name is None:
            import tempfile

            _, file_name = tempfile.mkstemp()

    n_jobs = effective_n_jobs(n_jobs)
    rep = update_rep(rep)

    if ("W_" + rep) not in data.uns:
        neighbors(
            data,
            K=K,
            rep=rep,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )

    indices_key = rep + "_knn_indices"
    distances_key = rep + "_knn_distances"

    if not knn_is_cached(data, indices_key, distances_key, select_K):
        raise ValueError("Please run neighbors first!")

    selected = select_cells(
        data.uns[distances_key],
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )

    X_full = X_from_rep(data, rep)
    X = X_full[selected, :]

    ds_indices_key = "ds_" + rep + "_knn_indices"
    ds_distances_key = "ds_" + rep + "_knn_distances"
    indices, distances = calculate_nearest_neighbors(
        X, K=K, n_jobs=n_jobs, random_state=random_state, full_speed=full_speed
    )
    data.uns[ds_indices_key] = indices
    data.uns[ds_distances_key] = distances

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

    n_components = 2 if not is3d else 3
    Y_init = np.zeros((data.shape[0], n_components), dtype=np.float64)
    Y_init[selected, :] = X_fle
    Y_init[~selected, :] = net_train_and_predict(
        X, X_fle, X_full[~selected, :], net_alpha, random_state, verbose=True
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    data.obsm["X_" + out_basis] = calc_force_directed_layout(
        W_from_rep(data, rep),
        file_name,
        n_jobs,
        target_change_per_node,
        polish_target_steps,
        is3d,
        memory,
        random_state,
        init=Y_init,
    )
