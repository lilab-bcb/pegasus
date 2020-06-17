import time
import numpy as np

from scipy.sparse import issparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from scipy.stats import entropy
from sklearn.decomposition import PCA
from sklearn.utils.extmath import randomized_svd
from typing import List, Tuple
from pegasusio import MultimodalData

from pegasus.tools import update_rep, W_from_rep

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


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


def calc_von_neumann_entropy(lambdas: List[float], t: float) -> float:
    etas = 1.0 - lambdas ** t
    etas = etas / etas.sum()
    return entropy(etas)


def find_knee_point(x: List[float], y: List[float]) -> int:
    """ Return the knee point, which is defined as the point furthest from line between two end points
    """
    p1 = np.array((x[0], y[0]))
    p2 = np.array((x[-1], y[-1]))
    length_p12 = np.linalg.norm(p2 - p1)

    max_dis = 0.0
    knee = 0
    for cand_knee in range(1, len(x) - 1):
        p3 = np.array((x[cand_knee], y[cand_knee]))
        dis = np.linalg.norm(np.cross(p2 - p1, p3 - p1)) / length_p12
        if max_dis < dis:
            max_dis = dis
            knee = cand_knee

    return knee


def calculate_diffusion_map(
    W: "csr_matrix", n_components: int, solver: str, random_state: int, max_t: int
) -> Tuple["np.array", "np.array", "np.array"]:
    assert issparse(W)

    nc, labels = connected_components(W, directed=True, connection="strong")
    logger.info("Calculating connected components is done.")

    assert nc == 1

    W_norm, diag, diag_half = calculate_normalized_affinity(W)
    logger.info("Calculating normalized affinity matrix is done.")

    if solver == "eigsh":
        np.random.seed(random_state)
        v0 = np.random.uniform(-1.0, 1.0, W_norm.shape[0])
        Lambda, U = eigsh(W_norm, k=n_components, v0=v0)
        Lambda = Lambda[::-1]
        U = U[:, ::-1]
    else:
        assert solver == "randomized"
        U, S, VT = randomized_svd(
            W_norm, n_components=n_components, random_state=random_state
        )
        signs = np.sign((U * VT.transpose()).sum(axis=0))  # get eigenvalue signs
        Lambda = signs * S  # get eigenvalues

    # remove the first eigen value and vector
    Lambda = Lambda[1:]
    U = U[:, 1:]
    Phi = U / diag_half[:, np.newaxis]

    if max_t == -1:
        Lambda_new = Lambda / (1.0 - Lambda)
    else:
        # Find the knee point
        x = np.array(range(1, max_t + 1), dtype = float)
        y = np.array([calc_von_neumann_entropy(Lambda, t) for t in x])
        t = x[find_knee_point(x, y)]
        logger.info("Detected knee point at t = {:.0f}.".format(t))

        # U_df = U * Lambda #symmetric diffusion component
        Lambda_new = Lambda * ((1.0 - Lambda ** t) / (1.0 - Lambda))
    Phi_pt = Phi * Lambda_new  # asym pseudo component

    return Phi_pt, Lambda, Phi  # , U_df, W_norm


@timer(logger=logger)
def diffmap(
    data: MultimodalData,
    n_components: int = 100,
    rep: str = "pca",
    solver: str = "eigsh",
    random_state: int = 0,
    max_t: float = 5000,
) -> None:
    """Calculate Diffusion Map.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``100``
        Number of diffusion components to calculate.

    rep: ``str``, optional, default: ``"pca"``
        Embedding Representation of data used for calculating the Diffusion Map. By default, use PCA coordinates.

    solver: ``str``, optional, default: ``"eigsh"``
        Solver for eigen decomposition:
            * ``"eigsh"``: default setting. Use *scipy* `eigsh <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html>`_ as the solver to find eigenvalus and eigenvectors using the Implicitly Restarted Lanczos Method.
            * ``"randomized"``: Use *scikit-learn* `randomized_svd <https://scikit-learn.org/stable/modules/generated/sklearn.utils.extmath.randomized_svd.html>`_ as the solver to calculate a truncated randomized SVD.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    max_t: ``float``, optional, default: ``5000``
        pegasus tries to determine the best t to sum up to between ``[1, max_t]``.

    Returns
    -------
    ``None``

    Update ``data.obsm``:
        * ``data.obsm["X_diffmap"]``: Diffusion Map matrix of the data.

    Update ``data.uns``:
        * ``data.uns["diffmap_evals"]``: Eigenvalues corresponding to Diffusion Map matrix.

    Examples
    --------
    >>> pg.diffmap(data)
    """

    rep = update_rep(rep)
    Phi_pt, Lambda, Phi = calculate_diffusion_map(
        W_from_rep(data, rep),
        n_components=n_components,
        solver=solver,
        random_state=random_state,
        max_t = max_t,
    )

    data.obsm["X_diffmap"] = Phi_pt
    data.uns["diffmap_evals"] = Lambda
    data.obsm["X_phi"] = Phi
    # data.uns['W_norm'] = W_norm
    # data.obsm['X_dmnorm'] = U_df


@timer(logger=logger)
def reduce_diffmap_to_3d(data: MultimodalData, random_state: int = 0) -> None:
    """Reduce high-dimensional Diffusion Map matrix to 3-dimentional.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
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
    >>> pg.reduce_diffmap_to_3d(data)
    """

    if "X_diffmap" not in data.obsm.keys():
        raise ValueError("Please run diffmap first!")

    pca = PCA(n_components=3, random_state=random_state)
    data.obsm["X_diffmap_pca"] = pca.fit_transform(data.obsm["X_diffmap"])
