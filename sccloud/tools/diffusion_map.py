import time
import numpy as np
import logging

from scipy.sparse import issparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from sklearn.decomposition import PCA
from sklearn.utils.extmath import randomized_svd
from typing import Tuple

from sccloud.tools import update_rep, W_from_rep

logger = logging.getLogger("sccloud")


def calculate_normalized_affinity(
    W: "csr_matrix"
) -> Tuple["csr_matrix", "np.array", "np.array"]:
    diag = W.sum(axis=1).A1
    diag_half = np.sqrt(diag)
    W_norm = W.tocoo(copy=True)
    W_norm.data /= diag_half[W_norm.row]
    W_norm.data /= diag_half[W_norm.col]
    W_norm = W_norm.tocsr()

    return W_norm, diag, diag_half


def calculate_diffusion_map(
    W: "csr_matrix", n_components: int, alpha: float, solver: str, random_state: int
) -> Tuple["np.array", "np.array"]:
    assert issparse(W)

    nc, labels = connected_components(W, directed=True, connection="strong")
    logger.info("Calculating connected components is done.")

    assert nc == 1

    W_norm, diag, diag_half = calculate_normalized_affinity(W)
    logger.info("Calculating normalized affinity matrix is done.")

    if solver == "randomized":
        U, S, VT = randomized_svd(
            W_norm, n_components=n_components, random_state=random_state
        )
        signs = np.sign((U * VT.transpose()).sum(axis=0))  # get eigenvalue signs
        Lambda = signs * S  # get eigenvalues
    else:
        assert solver == "eigsh"
        np.random.seed(random_state)
        v0 = np.random.uniform(-1.0, 1.0, W_norm.shape[0])
        Lambda, U = eigsh(W_norm, k=n_components, v0=v0)
        Lambda = Lambda[::-1]
        U = U[:, ::-1]

    # remove the first eigen value and vector
    Lambda = Lambda[1:]
    U = U[:, 1:]

    Phi = U / diag_half[:, np.newaxis]
    Lambda_new = Lambda / (1 - alpha * Lambda)

    # U_df = U * Lambda #symmetric diffusion component
    Phi_pt = Phi * Lambda_new  # asym pseudo component

    return Phi_pt, Lambda  # , U_df, W_norm


def diffmap(
    data: "AnnData",
    n_components: int = 50,
    rep: str = "pca",
    alpha: float = 0.5,
    solver: str = "randomized",
    random_state: int = 0,
) -> None:
    """Calculate Diffusion Map.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``50``
        Number of diffusion components to calculate.

    rep: ``str``, optional, default: ``"pca"``
        Embedding Representation of data used for calculating the Diffusion Map. By default, use PCA coordinates.

    alpha: ``float``, optional, default: ``0.5``
        Power parameter for diffusion-based pseudotime.

    solver: ``str``, optional, default: ``"randomized"``
        Solver for eigen decomposition:
            * ``"randomized"``: default setting. Use *scikit-learn* `randomized_svd <https://scikit-learn.org/stable/modules/generated/sklearn.utils.extmath.randomized_svd.html>`_ as the solver to calculate a truncated randomized SVD.
            * ``"eigsh"``: Use *scipy* `eigsh <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html>`_ as the solver to find eigenvalus and eigenvectors using the Implicitly Restarted Lanczos Method.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm["X_diffmap"]``: Diffusion Map matrix of the data.

    Update ``data.uns``:
        * ``data.uns["diffmap_evals"]``: Eigenvalues corresponding to Diffusion Map matrix.

    Examples
    --------
    >>> scc.diffmap(adata)
    """

    start = time.time()
    rep = update_rep(rep)
    Phi_pt, Lambda = calculate_diffusion_map(
        W_from_rep(data, rep),
        n_components=n_components,
        alpha=alpha,
        solver=solver,
        random_state=random_state,
    )

    data.obsm["X_diffmap"] = Phi_pt
    data.uns["diffmap_evals"] = Lambda
    # data.uns['W_norm'] = W_norm
    # data.obsm['X_dmnorm'] = U_df

    end = time.time()
    logger.info("diffmap finished. Time spent = {:.2f}s.".format(end - start))


def reduce_diffmap_to_3d(data: "AnnData", random_state: int = 0) -> None:
    """Reduce high-dimensional Diffusion Map matrix to 3-dimentional.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.


    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm["X_diffmap_pca"]``: 3D Diffusion Map matrix of data.

    Examples
    --------
    >>> scc.reduce_diffmap_to_3d(adata)
    """
    start = time.time()

    if "X_diffmap" not in data.obsm.keys():
        raise ValueError("Please run diffmap first!")

    pca = PCA(n_components=3, random_state=random_state)
    data.obsm["X_diffmap_pca"] = pca.fit_transform(data.obsm["X_diffmap"])

    end = time.time()
    logger.info(
        "Reduce diffmap to 3D is done. Time spent = {:.2f}s.".format(end - start)
    )
