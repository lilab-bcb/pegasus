import time
import numpy as np
import logging

from scipy.sparse import issparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from sklearn.decomposition import PCA
from sklearn.utils.extmath import randomized_svd
from typing import Tuple

from scCloud.tools import update_rep, W_from_rep

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
    """
    TODO: documentation.
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
    """
    TODO: documentation.
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
