import time
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype

from scipy.sparse import issparse
from collections import defaultdict
from joblib import Parallel, delayed
import skmisc.loess as sl
from typing import List
from pegasusio import MultimodalData

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


def estimate_feature_statistics(data: MultimodalData, consider_batch: bool) -> None:
    """ Estimate feature (gene) statistics per channel, such as mean, var etc.
    """
    assert issparse(data.X)

    if consider_batch:
        start = time.perf_counter()
        # The reason that we test if 'Channel' and 'Channels' exist in addition to the highly_variable_features function is for the case that we do not perform feature selection but do batch correction
        if "Channel" not in data.obs:
            data.obs["Channel"] = ""

        if "Channels" not in data.uns:
            data.uns["Channels"] = data.obs["Channel"].cat.categories.values if is_categorical_dtype(data.obs["Channel"]) else data.obs["Channel"].unique()

        if data.uns["Channels"].size == 1:
            return None

        if "Group" not in data.obs:
            data.obs["Group"] = "one_group"

        if "Groups" not in data.uns:
            data.uns["Groups"] = data.obs["Group"].cat.categories.values if is_categorical_dtype(data.obs["Group"]) else data.obs["Group"].unique()

        channels = data.uns["Channels"]
        groups = data.uns["Groups"]

        ncells = np.zeros(channels.size)
        means = np.zeros((data.shape[1], channels.size))
        partial_sum = np.zeros((data.shape[1], channels.size))

        group_dict = defaultdict(list)
        for i, channel in enumerate(channels):
            idx = np.isin(data.obs["Channel"], channel)
            mat = data.X[idx].astype(np.float64)
            ncells[i] = mat.shape[0]

            if ncells[i] == 0:
                continue

            if ncells[i] == 1:
                means[:, i] = mat.toarray()[0]
            else:
                means[:, i] = mat.mean(axis=0).A1
                m2 = mat.power(2).sum(axis=0).A1
                partial_sum[:, i] = m2 - ncells[i] * (means[:, i] ** 2)

            group = data.obs["Group"][idx.nonzero()[0][0]]
            group_dict[group].append(i)

        partial_sum[partial_sum < 1e-6] = 0.0

        overall_means = np.dot(means, ncells) / data.shape[0]
        batch_adjusted_vars = np.zeros(data.shape[1])

        c2gid = np.zeros(channels.size, dtype=int)
        gncells = np.zeros(groups.size)
        gmeans = np.zeros((data.shape[1], groups.size))
        gstds = np.zeros((data.shape[1], groups.size))

        for i, group in enumerate(groups):
            gchannels = group_dict[group]
            c2gid[gchannels] = i
            gncells[i] = ncells[gchannels].sum()
            gmeans[:, i] = np.dot(means[:, gchannels], ncells[gchannels]) / gncells[i]
            gstds[:, i] = (
                partial_sum[:, gchannels].sum(axis=1) / gncells[i]
            ) ** 0.5  # calculate std
            if groups.size > 1:
                batch_adjusted_vars += gncells[i] * (
                    (gmeans[:, i] - overall_means) ** 2
                )

        data.varm["means"] = means
        data.varm["partial_sum"] = partial_sum
        data.uns["ncells"] = ncells

        data.varm["gmeans"] = gmeans
        data.varm["gstds"] = gstds
        data.uns["gncells"] = gncells
        data.uns["c2gid"] = c2gid

        data.var["mean"] = overall_means
        data.var["var"] = (batch_adjusted_vars + partial_sum.sum(axis=1)) / (
            data.shape[0] - 1.0
        )
        end = time.perf_counter()
        logger.info(
            "Estimation on feature statistics per channel is finished. Time spent = {:.2f}s.".format(
                end - start
            )
        )
    else:
        mean = data.X.mean(axis=0).A1
        m2 = data.X.power(2).sum(axis=0).A1
        var = (m2 - data.X.shape[0] * (mean ** 2)) / (data.X.shape[0] - 1)

        data.var["mean"] = mean
        data.var["var"] = var


def fit_loess(x: List[float], y: List[float], span: float, degree: int) -> object:
    try:
        lobj = sl.loess(x, y, span=span, degree=2)
        lobj.fit()
        return lobj
    except ValueError:
        return None


def select_hvf_pegasus(
    data: MultimodalData, consider_batch: bool, n_top: int = 2000, span: float = 0.02
) -> None:
    """ Select highly variable features using the pegasus method
    """
    if "robust" not in data.var:
        raise ValueError("Please run `qc_metrics` to identify robust genes")

    estimate_feature_statistics(data, consider_batch)

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
    X: "csr_matrix",
    n_top: int,
    min_disp: float,
    max_disp: float,
    min_mean: float,
    max_mean: float,
) -> List[int]:
    """ HVF selection for one channel using Seurat method
    """
    X = X.copy().expm1()
    mean = X.mean(axis=0).A1
    m2 = X.power(2).sum(axis=0).A1
    var = (m2 - X.shape[0] * (mean ** 2)) / (X.shape[0] - 1)

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
    X: "csr_matrix",
    channels: List[str],
    cell2channel: List[str],
    n_top: int,
    n_jobs: int,
    min_disp: float,
    max_disp: float,
    min_mean: float,
    max_mean: float,
) -> List[int]:
    Xs = []
    for channel in channels:
        Xs.append(X[np.isin(cell2channel, channel)])

    from joblib import effective_n_jobs

    n_jobs = effective_n_jobs(n_jobs)

    res_arr = np.array(
        Parallel(n_jobs=n_jobs)(
            delayed(select_hvf_seurat_single)(
                Xs[i], n_top, min_disp, max_disp, min_mean, max_mean
            )
            for i in range(channels.size)
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
    data: MultimodalData,
    consider_batch: bool,
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
        select_hvf_seurat_multi(
            X,
            data.uns["Channels"],
            data.obs["Channel"],
            n_top,
            n_jobs=n_jobs,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
        )
        if consider_batch
        else select_hvf_seurat_single(
            X,
            n_top=n_top,
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
    data: MultimodalData,
    consider_batch: bool,
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

    consider_batch: ``bool``.
        Whether consider batch effects or not.

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
        Number of threads to be used during calculation. If ``-1``, all available threads will be used.

    Returns
    -------

    Examples
    --------
    >>> pg.highly_variable_features(data, consider_batch = False)
    """
    if "Channels" not in data.uns:
        if "Channel" not in data.obs:
            data.obs["Channel"] = ""
        data.uns["Channels"] = data.obs["Channel"].cat.categories.values if is_categorical_dtype(data.obs["Channel"]) else data.obs["Channel"].unique()

    if data.uns["Channels"].size == 1 and consider_batch:
        consider_batch = False
        logger.warning(
            "Warning: only contains one channel, no need to consider batch for selecting highly variable features."
        )

    if flavor == "pegasus":
        select_hvf_pegasus(data, consider_batch, n_top=n_top, span=span)
    else:
        assert flavor == "Seurat"
        select_hvf_seurat(
            data,
            consider_batch,
            n_top=n_top,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
            n_jobs=n_jobs,
        )

    logger.info(
        "{} highly variable features have been selected.".format(
            data.var["highly_variable_features"].sum()
        )
    )
