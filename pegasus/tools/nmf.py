import numpy as np
import pandas as pd

from typing import List, Tuple
from pegasusio import UnimodalData, MultimodalData
from pegasus.tools import eff_n_jobs, select_features

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer, run_gc



@timer(logger=logger)
def nmf(
    data: MultimodalData,
    n_components: int = 20,
    features: str = "highly_variable_features",
    max_value: float = None,
    space: str = "log",
    init: str = "nndsvdar",
    algo: str = "hals",
    mode: str = "batch",
    tol: float = 1e-4,
    use_gpu: bool = False,
    alpha_W: float = 0.0,
    l1_ratio_W: float = 0.0,
    alpha_H: float = 0.0,
    l1_ratio_H: float = 0.0,
    fp_precision: str = "float",
    n_jobs: int = -1,
    random_state: int = 0,
) -> None:
    """Perform Nonnegative Matrix Factorization (NMF) to the data using Frobenius norm.

    The calculation uses *NMF-Torch*.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``50``.
        Number of Principal Components to get.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for PCA.

    max_value: ``float``, optional, default: ``None``.
        The threshold to truncate data symmetrically after scaling. If ``None``, do not truncate.

    space: ``str``, optional, default: ``log``.
        Choose from ``log`` and ``expression``. ``log`` works on log-transformed expression space; ``expression`` works on the original expression space (normalized by total UMIs).

    init: ``str``, optional, default: ``nndsvdar``.
        Method to initialize NMF. Options are 'random', 'nndsvd', 'nndsvda' and 'nndsvdar'.

    algo: ``str``, optional, default: ``hals``
        Choose from ``mu`` (Multiplicative Update), ``hals`` (Hierarchical Alternative Least Square) and ``bpp`` (alternative non-negative least squares with Block Principal Pivoting method).

    mode: ``str``, optional, default: ``batch``
        Learning mode. Choose from ``batch`` and ``online``. Notice that ``online`` only works when ``beta=2.0``. For other beta loss, it switches back to ``batch`` method.        

    tol: ``float``, optional, default: ``1e-4``
        The toleration used for convergence check.

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

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``.
        Random seed to be set for reproducing result.

    Returns
    -------
    ``None``.

    Update ``data.obsm``:

        * ``data.obsm["X_nmf"]``: Scaled NMF coordinates. Each column has a unit variance.

    Update ``data.uns``:

        * ``data.uns["W"]``: The feature factor matrix. 

        * ``data.uns["H"]``: The coordinate factor matrix.

        * ``data.uns["nmf_features"]``: Record the features used to perform NMF analysis.

    Examples
    --------
    >>> pg.nmf(data)
    """
    keyword = select_features(data, features=features, scale=True, max_value=max_value, expression_space=space!="log")
    X = data.uns[keyword]
    X[X < 0] = 0.0

    try:
        from nmf import run_nmf
    except ImportError as e:
        import sys
        logger.error(f"{e}\nNeed Harmony! Try 'pip install harmony-pytorch'.")
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
    )

    data.uns["nmf_features"] = features # record which feature to use
    data.uns["W"] = np.ascontiguousarray(W.T, dtype=np.float32)  # cannot be varm because numbers of features are not the same
    data.uns["H"] = np.ascontiguousarray(H, dtype=np.float32)
    H = data.uns["H"]
    data.obsm["X_nmf"] = H / H.std(axis = 0) # scale to unit variance
