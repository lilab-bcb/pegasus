import time
import numpy as np
import scipy
import umap as umap_module
import uuid

from threadpoolctl import threadpool_limits
from pegasusio import MultimodalData
from pynndescent import NNDescent

from pegasus.tools import (
    eff_n_jobs,
    update_rep,
    X_from_rep,
    W_from_rep,
    get_neighbors,
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
    nthreads,
    no_dims,
    perplexity,
    early_exag_coeff,
    learning_rate,
    rand_seed,
    initialization=None,
    max_iter=750,
    stop_early_exag_iter=250,
    mom_switch_iter=250,
):
    """
    TODO: Calculate t-SNE embeddings using the FIt-SNE package
    """
    # FItSNE will change X content

    # Check if fftw3 is installed.
    import ctypes.util

    fftw3_loc = ctypes.util.find_library("fftw3")
    if fftw3_loc is None:
        raise Exception("Please install 'fftw3' first to use the FIt-SNE feature!")

    try:
        from fitsne import FItSNE
    except ModuleNotFoundError:
        print("Exception: FItSNE package is not installed, please install fitsne: pip install fitsne")
        exit()

    return FItSNE(
        X,
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

class DummyNNDescent(NNDescent):
    def __init__(self):
        None

# Running umap using our own kNN indices
def calc_umap(
    X,
    n_components,
    n_neighbors,
    min_dist,
    spread,
    random_state,
    densmap,
    dens_lambda,
    dens_frac,
    dens_var_shift,
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
        densmap=densmap,
        dens_lambda=dens_lambda,
        dens_frac=dens_frac,
        dens_var_shift=dens_var_shift,
        verbose=True,
    )

    if X.shape[0] < 4096 or knn_indices is None:
        logger.info(f"Using umap kNN graph because number of cells {X.shape[0]} is smaller than 4096 or knn_indices is not provided.")
    else:
        assert knn_dists is not None, "No kNN graph is found! Please calculate it by 'pegasus.neighbors' function first!"
        dummy_nnd = DummyNNDescent()
        umap_obj.precomputed_knn = (knn_indices, knn_dists, dummy_nnd)

    return umap_obj.fit_transform(X)


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
    try:
        import forceatlas2 as fa2
    except ModuleNotFoundError:
        print("Exception: Please install module forceatlas2-python: pip install forceatlas2-python")
        exit()
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
    learning_rate: float = "auto",
    initialization: str = "pca",
    random_state: int = 0,
    out_basis: str = "tsne",
) -> None:
    """Calculate t-SNE embedding of cells using the FIt-SNE package.

    This function uses fitsne_ package. See [Linderman19]_ for details on FIt-SNE algorithm.

    .. _fitsne: https://github.com/KlugerLab/FIt-SNE

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated FI-tSNE coordinates. By default, generate 2-dimensional data for 2D visualization.

    perplexity: ``float``, optional, default: ``30``
        The perplexity is related to the number of nearest neighbors used in other manifold learning algorithms. Larger datasets usually require a larger perplexity.

    early_exaggeration: ``int``, optional, default: ``12``
        Controls how tight natural clusters in the original space are in the embedded space, and how much space will be between them.

    learning_rate: ``float``, optional, default: ``auto``
        By default, the learning rate is determined automatically as max(data.shape[0] / early_exaggeration, 200). See [Belkina19]_ and [Kobak19]_ for details.

    initialization: ``str``, optional, default: ``pca``
        Initialization can be either ``pca`` or ``random`` or np.ndarray. By default, we use ``pca`` initialization according to [Kobak19]_.

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
    >>> pg.tsne(data)
    """

    rep = update_rep(rep)
    n_jobs = eff_n_jobs(n_jobs)
    X = X_from_rep(data, rep).astype(np.float64)

    if learning_rate == "auto":
        learning_rate = max(X.shape[0] / early_exaggeration, 200.0)

    if initialization == "random":
        initialization = None
    elif initialization == "pca":
        if rep == "pca":
            initialization = X[:, 0:n_components].copy()
        else:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=n_components, random_state=random_state)
            with threadpool_limits(limits = n_jobs):
                initialization = np.ascontiguousarray(pca.fit_transform(X))
        initialization = initialization / np.std(initialization[:, 0]) * 0.0001
    else:
        assert isinstance(initialization, np.ndarray) and initialization.ndim == 2 and initialization.shape[0] == X.shape[0] and initialization.shape[1] == n_components
        if initialization.dtype != np.float64:
            initialization = initialization.astype(np.float64)

    key = f"X_{out_basis}"
    data.obsm[key] = calc_tsne(
        X,
        n_jobs,
        n_components,
        perplexity,
        early_exaggeration,
        learning_rate,
        random_state,
        initialization,
    )
    data.register_attr(key, "basis")


@timer(logger=logger)
def umap(
    data: MultimodalData,
    rep: str = "pca",
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.5,
    spread: float = 1.0,
    densmap: bool = False,
    dens_lambda: float = 2.0,
    dens_frac: float = 0.3,
    dens_var_shift: float = 0.1,
    n_jobs: int = -1,
    full_speed: bool = False,
    random_state: int = 0,
    out_basis: str = "umap",
) -> None:
    """Calculate UMAP embedding of cells.

    This function uses umap-learn_ package. See [McInnes18]_ for details on UMAP.

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

    densmap: ``bool``, optional, default: ``False``
        Whether the density-augmented objective of densMAP should be used for optimization, which will generate an embedding where
        local densities are encouraged to be correlated with those in the original space.

    dens_lambda: ``float``, optional, default: ``2.0``
        Controls the regularization weight of the density correlation term in densMAP. Only works when *densmap* is ``True``.
        Larger values prioritize density preservation over the UMAP objective, while values closer to 0 for the opposite direction.
        Notice that setting this parameter to ``0`` is equivalent to running the original UMAP algorithm.

    dens_frac: ``float``, optional, default: ``0.3``
        Controls the fraction of epochs (between 0 and 1) where the density-augmented objective is used in densMAP. Only works when
        *densmap* is ``True``.
        The first ``(1 - dens_frac)`` fraction of epochs optimize the original UMAP objective before introducing the density
        correlation term.

    dens_var_shift: ``float``, optional, default, ``0.1``
        A small constant added to the variance of local radii in the embedding when calculating the density correlation objective to
        prevent numerical instability from dividing by a small number. Only works when *densmap* is ``True``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use for computing kNN graphs. If ``-1``, use all physical CPU cores.

    full_speed: ``bool``, optional, default: ``False``
        * If ``True``, use multiple threads in constructing ``hnsw`` index. However, the kNN results are not reproducible.
        * Otherwise, use only one thread to make sure results are reproducible.

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
    rep = update_rep(rep)
    X = X_from_rep(data, rep)

    if data.shape[0] < n_neighbors:
        logger.warning(f"Warning: Number of samples = {data.shape[0]} < K = {n_neighbors}!\n Set K to {data.shape[0]}.")
        n_neighbors = data.shape[0]

    knn_indices, knn_dists = get_neighbors(data, K = n_neighbors, rep = rep, n_jobs = n_jobs, random_state = random_state, full_speed = full_speed)
    knn_indices = np.insert(knn_indices[:, 0 : n_neighbors - 1], 0, range(data.shape[0]), axis=1)
    knn_dists = np.insert(knn_dists[:, 0 : n_neighbors - 1], 0, 0.0, axis=1)

    key = f"X_{out_basis}"
    data.obsm[key] = calc_umap(
        X,
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        densmap=densmap,
        dens_lambda=dens_lambda,
        dens_frac=dens_frac,
        dens_var_shift=dens_var_shift,
        random_state=random_state,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )
    data.register_attr(key, "basis")


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
    """Construct the Force-directed (FLE) graph.

    This implementation uses forceatlas2-python_ package, which is a Python wrapper of ForceAtlas2_.

    See [Jacomy14]_ for details on FLE.

    .. _forceatlas2-python: https://github.com/klarman-cell-observatory/forceatlas2-python
    .. _ForceAtlas2: https://github.com/klarman-cell-observatory/forceatlas2

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    file_name: ``str``, optional, default: ``None``
        Temporary file to store the coordinates as the input to forceatlas2. If ``None``, use ``tempfile.mkstemp`` to generate file name.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

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

    rep = update_rep(rep)
    n_jobs = eff_n_jobs(n_jobs)

    if ("W_" + rep) not in data.uns:
        neighbors(
            data,
            K=K,
            rep=rep,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )

    key = f"X_{out_basis}"
    data.obsm[key] = calc_force_directed_layout(
        W_from_rep(data, rep),
        file_name,
        n_jobs,
        target_change_per_node,
        target_steps,
        is3d,
        memory,
        random_state,
    )
    data.register_attr(key, "basis")


@timer(logger=logger)
def select_cells(distances, frac, K=25, alpha=1.0, random_state=0):
    """
    TODO: documentation (not user API)
    """
    nsample = distances.shape[0]
    assert K >= 2
    if K > distances.shape[1] + 1:
        logger.info(f"Warning: in select_cells, K = {K} > the number of calculated nearest neighbors {distances.shape[1] + 1}!\nSet K to {distances.shape[1] + 1}")
        K = distances.shape[1] + 1

    probs = np.zeros(nsample)
    if alpha == 0.0:
        probs[:] = 1.0  # uniform
    elif alpha == 1.0:
        probs[:] = distances[:, K - 2]
    else:
        probs[:] = distances[:, K - 2] ** alpha
    probs /= probs.sum()

    np.random.seed(random_state)
    selected = np.zeros(nsample, dtype=bool)
    selected[
        np.random.choice(nsample, size=int(nsample * frac), replace=False, p=probs)
    ] = True

    return selected


@timer(logger=logger)
def net_umap(
    data: MultimodalData,
    rep: str = "pca",
    n_jobs: int = -1,
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.5,
    spread: float = 1.0,
    densmap: bool = False,
    dens_lambda: float = 2.0,
    dens_frac: float = 0.3,
    dens_var_shift: float = 0.1,
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
    """Calculate Net-UMAP embedding of cells.

    Net-UMAP is an approximated UMAP embedding using Deep Learning model to improve the speed.

    In specific, the deep model used is MLPRegressor_, the *scikit-learn* implementation of Multi-layer Perceptron regressor.

    See [Li20]_ for details.

    .. _MLPRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    rep: ``str``, optional, default: ``"pca"``
        Representation of data used for the calculation. By default, use PCA coordinates. If ``None``, use the count matrix ``data.X``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

    n_components: ``int``, optional, default: ``2``
        Dimension of calculated UMAP coordinates. By default, generate 2-dimensional data for 2D visualization.

    n_neighbors: ``int``, optional, default: ``15``
        Number of nearest neighbors considered during the computation.

    min_dist: ``float``, optional, default: ``0.5``
        The effective minimum distance between embedded data points.

    spread: ``float``, optional, default: ``1.0``
        The effective scale of embedded data points.

    densmap: ``bool``, optional, default: ``False``
        Whether the density-augmented objective of densMAP should be used for optimization, which will generate an embedding where
        local densities are encouraged to be correlated with those in the original space.

    dens_lambda: ``float``, optional, default: ``2.0``
        Controls the regularization weight of the density correlation term in densMAP. Only works when *densmap* is ``True``.
        Larger values prioritize density preservation over the UMAP objective, while values closer to 0 for the opposite direction.
        Notice that setting this parameter to ``0`` is equivalent to running the original UMAP algorithm.

    dens_frac: ``float``, optional, default: ``0.3``
        Controls the fraction of epochs (between 0 and 1) where the density-augmented objective is used in densMAP. Only works when
        *densmap* is ``True``.
        The first ``(1 - dens_frac)`` fraction of epochs optimize the original UMAP objective before introducing the density
        correlation term.

    dens_var_shift: ``float``, optional, default, ``0.1``
        A small constant added to the variance of local radii in the embedding when calculating the density correlation objective to
        prevent numerical instability from dividing by a small number. Only works when *densmap* is ``True``.

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
    n_jobs = eff_n_jobs(n_jobs)
    knn_indices, knn_dists = get_neighbors(data, K = select_K, rep = rep, n_jobs = n_jobs, random_state = random_state, full_speed = full_speed)

    selected = select_cells(
        knn_dists,
        select_frac,
        K=select_K,
        alpha=select_alpha,
        random_state=random_state,
    )
    X_full = X_from_rep(data, rep)
    X = X_full[selected, :]

    if data.shape[0] < n_neighbors:
        logger.warning(f"Warning: Number of samples = {data.shape[0]} < K = {n_neighbors}!\n Set K to {data.shape[0]}.")
        n_neighbors = data.shape[0]

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
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        densmap=densmap,
        dens_lambda=dens_lambda,
        dens_frac=dens_frac,
        dens_var_shift=dens_var_shift,
        random_state=random_state,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )

    data.uns["X_" + out_basis + "_small"] = X_umap
    data.obs["ds_selected"] = selected

    Y_init = np.zeros((data.shape[0], n_components), dtype=np.float64)
    Y_init[selected, :] = X_umap
    Y_init[~selected, :] = net_train_and_predict(
        X, X_umap, X_full[~selected, :], net_alpha, n_jobs, random_state, verbose=True
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    knn_indices, knn_dists = get_neighbors(data, K = n_neighbors, rep = rep, n_jobs = n_jobs, random_state = random_state, full_speed = full_speed)
    knn_indices = np.insert(knn_indices[:, 0 : n_neighbors - 1], 0, range(data.shape[0]), axis=1)
    knn_dists = np.insert(knn_dists[:, 0 : n_neighbors - 1], 0, 0.0, axis=1)

    key = f"X_{out_basis}"
    data.obsm[key] = calc_umap(
        X_full,
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        spread=spread,
        densmap=densmap,
        dens_lambda=dens_lambda,
        dens_frac=dens_frac,
        dens_var_shift=dens_var_shift,
        random_state=random_state,
        init=Y_init,
        n_epochs=polish_n_epochs,
        learning_rate=polish_learning_rate,
        knn_indices=knn_indices,
        knn_dists=knn_dists,
    )
    data.register_attr(key, "basis")


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
    """Construct Net-Force-directed (FLE) graph.

    Net-FLE is an approximated FLE graph using Deep Learning model to improve the speed.

    In specific, the deep model used is MLPRegressor_, the *scikit-learn* implementation of Multi-layer Perceptron regressor.

    See [Li20]_ for details.

    .. _MLPRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPRegressor.html

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    file_name: ``str``, optional, default: ``None``
        Temporary file to store the coordinates as the input to forceatlas2. If ``None``, use ``tempfile.mkstemp`` to generate file name.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

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

    rep = update_rep(rep)
    n_jobs = eff_n_jobs(n_jobs)

    if ("W_" + rep) not in data.uns:
        neighbors(
            data,
            K=K,
            rep=rep,
            n_jobs=n_jobs,
            random_state=random_state,
            full_speed=full_speed,
        )

    knn_indices, knn_dists = get_neighbors(data, K = select_K, rep = rep, n_jobs = n_jobs, random_state = random_state, full_speed = full_speed)

    selected = select_cells(
        knn_dists,
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
        X, X_fle, X_full[~selected, :], net_alpha, n_jobs, random_state, verbose=True
    )

    data.obsm["X_" + out_basis + "_pred"] = Y_init

    key = f"X_{out_basis}"
    data.obsm[key] = calc_force_directed_layout(
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
    data.register_attr(key, "basis")
