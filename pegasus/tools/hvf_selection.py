import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from collections import defaultdict
from joblib import Parallel, delayed, parallel_backend
import skmisc.loess as sl
from typing import List, Union
from pegasusio import UnimodalData, MultimodalData

from pegasus.tools import eff_n_jobs, calc_mean_and_var, calc_expm1, calc_stat_per_batch, check_batch_key

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



@timer(logger=logger)
def estimate_feature_statistics(data: Union[MultimodalData, UnimodalData], batch: str) -> None:
    """ Estimate feature (gene) statistics per channel, such as mean, var etc.
    """
    if batch is None:
        data.var["mean"], data.var["var"] = calc_mean_and_var(data.X, axis=0)
    else:
        ncells, means, partial_sum = calc_stat_per_batch(data.X, data.obs[batch].values)
        partial_sum[partial_sum < 1e-6] = 0.0

        data.uns["ncells"] = ncells
        data.varm["means"] = means
        data.varm["partial_sum"] = partial_sum

        data.var["mean"] = np.dot(means, ncells) / data.shape[0]
        data.var["var"] = partial_sum.sum(axis=1) / (data.shape[0] - 1.0)
        

def fit_loess(x: List[float], y: List[float], span: float, degree: int) -> object:
    try:
        lobj = sl.loess(x, y, span=span, degree=degree)
        lobj.fit()
        return lobj
    except ValueError:
        return None


def select_hvf_pegasus(
    data: Union[MultimodalData, UnimodalData], batch: str, n_top: int = 2000, span: float = 0.02
) -> None:
    """ Select highly variable features using the pegasus method
    """
    if "robust" not in data.var:
        raise ValueError("Please run `identify_robust_genes` to identify robust genes")

    estimate_feature_statistics(data, batch)

    robust_idx = data.var["robust"].values
    hvf_index = np.zeros(robust_idx.sum(), dtype=bool)

    mean = data.var.loc[robust_idx, "mean"]
    var = data.var.loc[robust_idx, "var"]

    span_value = span
    while True:
        lobj = fit_loess(mean, var, span = span_value, degree = 2)
        if lobj is not None:
            break
        span_value += 0.01
    if span_value > span:
        logger.warning("Leoss span is adjusted from {:.2f} to {:.2f} to avoid fitting errors.".format(span, span_value))

    rank1 = np.zeros(hvf_index.size, dtype=int)
    rank2 = np.zeros(hvf_index.size, dtype=int)

    delta = var - lobj.outputs.fitted_values
    fc = var / lobj.outputs.fitted_values

    rank1[np.argsort(delta)[::-1]] = range(hvf_index.size)
    rank2[np.argsort(fc)[::-1]] = range(hvf_index.size)
    hvf_rank = rank1 + rank2

    hvf_index[np.argsort(hvf_rank)[:n_top]] = True

    data.var["hvf_loess"] = 0.0
    data.var.loc[robust_idx, "hvf_loess"] = lobj.outputs.fitted_values

    data.var["hvf_rank"] = -1
    data.var.loc[robust_idx, "hvf_rank"] = hvf_rank
    data.var["highly_variable_features"] = False
    data.var.loc[robust_idx, "highly_variable_features"] = hvf_index


def select_hvf_seurat_single(
    X: Union[csr_matrix, np.ndarray],
    n_top: int,
    min_disp: float,
    max_disp: float,
    min_mean: float,
    max_mean: float,
) -> List[int]:
    """ HVF selection for one channel using Seurat method
    """
    X = calc_expm1(X)

    mean, var = calc_mean_and_var(X, axis=0)

    dispersion = np.full(X.shape[1], np.nan)
    idx_valid = (mean > 0.0) & (var > 0.0)
    dispersion[idx_valid] = var[idx_valid] / mean[idx_valid]

    mean = np.log1p(mean)
    dispersion = np.log(dispersion)

    df = pd.DataFrame({"log_dispersion": dispersion, "bin": pd.cut(mean, bins=20)})
    log_disp_groups = df.groupby("bin")["log_dispersion"]
    log_disp_mean = log_disp_groups.mean()
    log_disp_std = log_disp_groups.std(ddof=1)
    log_disp_zscore = (
        df["log_dispersion"].values - log_disp_mean.loc[df["bin"]].values
    ) / log_disp_std.loc[df["bin"]].values
    log_disp_zscore[np.isnan(log_disp_zscore)] = 0.0

    hvf_rank = np.full(X.shape[1], -1, dtype=int)
    ords = np.argsort(log_disp_zscore)[::-1]

    if n_top is None:
        hvf_rank[ords] = range(X.shape[1])
        idx = np.logical_and.reduce(
            (
                mean > min_mean,
                mean < max_mean,
                log_disp_zscore > min_disp,
                log_disp_zscore < max_disp,
            )
        )
        hvf_rank[~idx] = -1
    else:
        hvf_rank[ords[:n_top]] = range(n_top)

    return hvf_rank


def select_hvf_seurat_multi(
    X: Union[csr_matrix, np.ndarray],
    batches: List[str],
    cell2batch: List[str],
    n_top: int,
    n_jobs: int,
    min_disp: float,
    max_disp: float,
    min_mean: float,
    max_mean: float,
) -> List[int]:
    Xs = []
    for batch in batches:
        Xs.append(X[np.isin(cell2batch, batch)])

    n_jobs = eff_n_jobs(n_jobs)
    with parallel_backend("loky", inner_max_num_threads=1):
        res_arr = np.array(
            Parallel(n_jobs=n_jobs)(
                delayed(select_hvf_seurat_single)(
                    Xs[i], n_top, min_disp, max_disp, min_mean, max_mean
                )
                for i in range(batches.size)
            )
        )

    selected = res_arr >= 0
    shared = selected.sum(axis=0)
    cands = (shared > 0).nonzero()[0]
    import numpy.ma as ma

    median_rank = ma.median(ma.masked_array(res_arr, mask=~selected), axis=0).data
    cands = sorted(cands, key=lambda x: median_rank[x])
    cands = sorted(cands, key=lambda x: shared[x], reverse=True)

    hvf_rank = np.full(X.shape[1], -1, dtype=int)
    hvf_rank[cands[:n_top]] = range(n_top)

    return hvf_rank


def select_hvf_seurat(
    data: Union[MultimodalData, UnimodalData],
    batch: str,
    n_top: int,
    min_disp: float,
    max_disp: float,
    min_mean: float,
    max_mean: float,
    n_jobs: int,
) -> None:
    """ Select highly variable features using Seurat method.
    """

    robust_idx = data.var["robust"].values
    X = data.X[:, robust_idx]

    hvf_rank = (
        select_hvf_seurat_single(
            X,
            n_top=n_top,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )
        if batch is None
        else select_hvf_seurat_multi(
            X,
            data.obs[batch].cat.categories.values,
            data.obs[batch],
            n_top,
            n_jobs=n_jobs,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )
    )

    hvf_index = hvf_rank >= 0

    data.var["hvf_rank"] = -1
    data.var.loc[robust_idx, "hvf_rank"] = hvf_rank
    data.var["highly_variable_features"] = False
    data.var.loc[robust_idx, "highly_variable_features"] = hvf_index


@timer(logger=logger)
def highly_variable_features(
    data: Union[MultimodalData, UnimodalData],
    batch: str = None,
    flavor: str = "pegasus",
    n_top: int = 2000,
    span: float = 0.02,
    min_disp: float = 0.5,
    max_disp: float = np.inf,
    min_mean: float = 0.0125,
    max_mean: float = 7,
    n_jobs: int = -1,
) -> None:
    """ Highly variable features (HVF) selection. The input data should be logarithmized.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    batch: ``str``, optional, default: ``None``
        A key in data.obs specifying batch information. If `batch` is not set, do not consider batch effects in selecting highly variable features. Otherwise, if `data.obs[batch]` is not categorical, `data.obs[batch]` will be automatically converted into categorical before highly variable feature selection.

    flavor: ``str``, optional, default: ``"pegasus"``
        The HVF selection method to use. Available choices are ``"pegasus"`` or ``"Seurat"``.

    n_top: ``int``, optional, default: ``2000``
        Number of genes to be selected as HVF. if ``None``, no gene will be selected.

    span: ``float``, optional, default: ``0.02``
        Only applicable when ``flavor`` is ``"pegasus"``. The smoothing factor used by *scikit-learn loess* model in pegasus HVF selection method.

    min_disp: ``float``, optional, default: ``0.5``
        Minimum normalized dispersion.

    max_disp: ``float``, optional, default: ``np.inf``
        Maximum normalized dispersion. Set it to ``np.inf`` for infinity bound.

    min_mean: ``float``, optional, default: ``0.0125``
        Minimum mean.

    max_mean: ``float``, optional, default: ``7``
        Maximum mean.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to be used during calculation. If ``-1``, all physical CPU cores will be used.

    Returns
    -------
    ``None``

    Update ``data.var``:
        * ``highly_variable_features``: replace with Boolean type array indicating the selected highly variable features.

    Examples
    --------
    >>> pg.highly_variable_features(data, consider_batch = False)
    """
    check_batch_key(data, batch, "Switch to not not consider batch effects for selecting highly variable features.")

    if flavor == "pegasus":
        select_hvf_pegasus(data, batch, n_top=n_top, span=span)
    else:
        assert flavor == "Seurat"
        select_hvf_seurat(
            data,
            batch,
            n_top=n_top,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
            n_jobs=n_jobs,
        )

    data.uns.pop("_tmp_fmat_highly_variable_features", None) # Pop up cached feature matrix

    logger.info(f"{data.var['highly_variable_features'].sum()} highly variable features have been selected.")
