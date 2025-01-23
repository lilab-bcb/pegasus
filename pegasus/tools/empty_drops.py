import itertools
import numba

import numpy as np
import pandas as pd
import scipy.optimize as so

from pegasusio import MultimodalData, UnimodalData
from pegasus.tools import eff_n_jobs
from pegasus.tools import SimpleGoodTuring
from pegasus.cylib.fast_utils import test_empty_drops
from pegasus.plotting import plot_barcode_rank

from anndata import AnnData
from scipy.special import loggamma
from scipy.sparse import csr_matrix, issparse
from statsmodels.stats.multitest import fdrcorrection as fdr
from joblib import Parallel, delayed, parallel_backend
from typing import Union, Optional

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


@timer(logger=logger)
def empty_drops(
    data: Union[MultimodalData, UnimodalData, AnnData],
    mat_key: Optional[str] = None,
    thresh_low: float = 100,
    distribution: str = "multi-dirichlet",
    alpha: Optional[float] = None,
    n_iters: int = 10000,
    random_state: int = 0,
    exclude_from: int = 50,
    significance_level: float = 0.01,
    show_summary_plot: bool = False,
    knee_method: str = "linear",
    n_jobs: int = 1,
) -> (int, int):
    if mat_key is None:
        if "counts" in data._unidata.matrices:
            mat_key = "counts"
        elif "raw.X" in data._unidata.matrices:
            mat_key = "raw.X"
    X = data.get_matrix(mat_key) if mat_key else data.X

    if "n_counts" not in data.obs:
        logger.info(
            "Calculate n_counts for EmptyDrops since 'n_counts' not in data.obs."
        )
        n_counts = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
    else:
        n_counts = data.obs["n_counts"].values

    assert (
        np.sum(X.sum(axis=0) == 0) == 0
    ), "Some genes have zero counts in data! Please run pegasus.identify_robust_genes() first!"

    idx_low = np.where(n_counts <= thresh_low)[0]
    G = X[idx_low, :]

    ambient_profile = G.sum(axis=0).A1
    sgt = SimpleGoodTuring(ambient_profile)
    ambient_proportion = sgt.get_proportions()

    use_alpha = True if distribution == "multi-dirichlet" else False

    if use_alpha:
        if alpha is None:
            alpha = _estimate_alpha(G, ambient_proportion)
            logger.info(f"Calculating alpha is finished. Estimate alpha = {alpha}.")
        else:
            logger.info(f"Use user-specified alpha = {alpha}.")

    idx_test = np.where(n_counts > thresh_low)[0]
    cells_test = data.obs.iloc[idx_test].index.values
    logger.info(f"Run test through {idx_test.size} cells with n_counts > {thresh_low}.")

    T = X[idx_test, :]
    tb = n_counts[idx_test]
    Lb_data = _logL_data_dep(T, ambient_proportion, alpha, use_alpha)
    Lb_alpha = _logL_data_indep(tb, alpha, use_alpha)

    n_jobs = eff_n_jobs(n_jobs)
    n_below = _test_empty_drops(use_alpha, alpha, ambient_proportion, tb, Lb_data, n_iters, random_state, n_jobs)
    pval = (n_below + 1) / (n_iters + 1)
    logger.info("Calculation on p-values is finished.")

    assert knee_method in ["spline", "linear"], "knee_method must be chosen from ['spline', 'linear']!"
    df_barcode_rank, knee, inflection = _rank_barcode(n_counts, thresh_low, exclude_from, knee_method)

    idx_always = np.where(tb >= knee)[0]
    pval_modified = pval.copy()
    pval_modified[idx_always] = 0.0
    logger.info(f"Adjust p-values by the knee point {knee}.")

    _, qval = fdr(pval_modified)

    data.obs["EmptyDrops.pval"] = np.nan
    data.obs.loc[cells_test, "EmptyDrops.pval"] = pval
    data.obs["EmptyDrops.qval"] = np.nan
    data.obs.loc[cells_test, "EmptyDrops.qval"] = qval
    data.obs["EmptyDrops.logL"] = np.nan
    data.obs.loc[cells_test, "EmptyDrops.logL"] = Lb_alpha + Lb_data
    data.obs["ambient"] = True
    data.obs.loc[data.obs["EmptyDrops.qval"] <= significance_level, "ambient"] = False

    if show_summary_plot:
        idx_nonambient = data.obs_names.get_indexer(data.obs.loc[~data.obs["ambient"]].index)
        plot_barcode_rank(df_barcode_rank, n_counts[idx_nonambient], thresh_low, knee, inflection)

    return (df_barcode_rank, knee, inflection)


def _logL_data_dep(Y, prop, alpha, use_alpha):
    if issparse(Y):
        idx_b, idx_g = Y.nonzero()
        y = Y.data
    else:
        idx_b, idx_g = np.nonzero(Y)
        y = Y[idx_b, idx_g]

    if use_alpha:
        alpha_prop = alpha * prop[idx_g]
        L_a1 = loggamma(y + alpha_prop) - loggamma(y + 1) - loggamma(alpha_prop)
    else:
        L_a1 = y * np.log(prop[idx_g]) - loggamma(y + 1)

    L_mat = csr_matrix((L_a1, (idx_b, idx_g)), shape=Y.shape)
    return L_mat.sum(axis=1).A1


def _logL_data_indep(tb, alpha, use_alpha):
    if use_alpha:
        return loggamma(tb + 1) + loggamma(alpha) - loggamma(tb + alpha)
    else:
        return loggamma(tb + 1)


def _rank_barcode(n_counts, thresh_low, exclude_from, knee_method="spline"):
    rank_values, rank_sizes = np.unique(n_counts, return_counts=True)
    rank_values = rank_values[::-1]
    rank_sizes = rank_sizes[::-1]
    tie_rank = np.cumsum(rank_sizes) - (rank_sizes - 1) / 2  # Get mid-rank of each tie

    idx_keep = np.where(rank_values > thresh_low)[0]
    y = np.log(rank_values[idx_keep])
    x = np.log(tie_rank[idx_keep])
    left_edge, right_edge = _find_curve_bounds(x, y, exclude_from)

    inflection = int(np.exp(y[right_edge]))

    fitted = np.full(tie_rank.size, np.nan)
    idx_focus = np.arange(left_edge, right_edge+1)
    if idx_focus.size >= 4:
        knee = _find_knee_point(x, y, fitted, idx_focus, method=knee_method)
        print(f"knee = {knee}")
    else:
        knee = int(np.ceil(np.exp(y[idx_focus[0]])))

    def repeat(vals, lens):
        return [element for element, count in zip(vals, lens) for _ in range(count)]

    df_rank = pd.DataFrame({
        "rank": tie_rank,
        "size": rank_sizes,
        "n_counts_obs": rank_values,
        "n_counts_fitted": fitted,
    })
    return (df_rank, knee, inflection)


def _find_knee_point(x, y, fitted, idx_focus, method="spline"):
    x_obs = x[idx_focus]
    y_obs = y[idx_focus]
    if method == "spline":
        from scipy.interpolate import UnivariateSpline

        spline = UnivariateSpline(x_obs, y_obs, k=3, s=5)
        fitted[idx_focus] = spline(x_obs)
        d1 = spline.derivative(n=1)(x_obs)
        d2 = spline.derivative(n=2)(x_obs)
        curvature = d2 / (1 + d1**2)**1.5
        knee = int(np.ceil(np.exp(y_obs[np.argmin(curvature)])))
    else:
        slope = (y_obs[-1] - y_obs[0]) / (x_obs[-1] - x_obs[0])
        intercept = y_obs[0] - x_obs[0] * slope
        y_fitted = x_obs * slope + intercept
        fitted[idx_focus] = y_fitted
        above = np.where(y_obs >= y_fitted)[0]
        distance = (y_obs[above] - y_fitted[above]) / np.sqrt(slope**2 + 1)
        knee = int(np.ceil(np.exp(y_obs[above[np.argmax(distance)]])))

    return knee


def _find_curve_bounds(x, y, exclude_from):
    d1n = np.diff(y) / np.diff(x)

    skip = np.min([d1n.size - 1, np.sum(x<=np.log(exclude_from))])
    d1n = d1n[-(d1n.size - skip):]

    right_edge = np.argmin(d1n)   # point with the least steep
    left_edge = np.argmax(d1n[:(right_edge+1)])   # point with the highest steep to the left of right edge

    return (left_edge + skip, right_edge + skip)


def _test_empty_drops(use_alpha, alpha, prop, tb, P_data, n_iters, random_state, n_jobs, temp_folder=None):
    idx_sorted = np.lexsort((P_data, tb))
    P_sorted = P_data[idx_sorted]
    tb_unique, tb_cnt = np.unique(tb[idx_sorted], return_counts=True)

    alpha_prop = alpha * prop if use_alpha else prop
    n_cells = tb.size
    n_genes = prop.size
    tb_max = tb_unique[-1]

    chunk_size = n_iters // n_jobs
    remainder = n_iters % n_jobs
    if chunk_size == 0:
        n_jobs = 1
        chunk_size = n_iters
        remainder = 0
    intervals = []
    start_pos = end_pos = 0
    for i in range(n_jobs):
        end_pos = start_pos + chunk_size + (i < remainder)
        if end_pos == start_pos:
            break
        intervals.append((start_pos, end_pos))
        start_pos = end_pos

    with parallel_backend("loky", inner_max_num_threads=1):
        result = Parallel(n_jobs=n_jobs, temp_folder=temp_folder)(
            delayed(test_empty_drops)(
                use_alpha,
                alpha_prop,
                tb_unique,
                tb_unique.size,
                tb_max,
                tb_cnt,
                P_sorted,
                n_cells,
                n_genes,
                random_state,
                intervals[i][0],
                intervals[i][1],
            )
            for i in range(n_jobs)
        )
    logger.info("Significance test is finished.")
    n_below_sorted = np.vstack(result).sum(axis=0)

    idx_inv = np.empty_like(idx_sorted)
    idx_inv[idx_sorted] = np.arange(idx_sorted.size)
    n_below = n_below_sorted[idx_inv]

    return n_below


def _calc_logL_by_increment(L1, alpha, prop, t, n_extra_counts, z, idx_extra_genes):
    return (
        L1
        + n_extra_counts * np.log(t / (t + alpha - 1))
        + np.sum(np.log((z[idx_extra_genes] + alpha * prop[idx_extra_genes] - 1) / z[idx_extra_genes]))
    )


def _estimate_alpha(G, prop, bounds=(0.01, 10000)):
    if issparse(G):
        tb = G.sum(axis=1).A1
        idx_g = G.nonzero()[1]
        y = G.data
    else:
        tb = G.sum(axis=1)
        idx_b, idx_g = np.nonzero(G)[1]
        y = G[idx_b, idx_g]

    n_cells = tb.size

    def neg_logL(alpha):
        # Remove terms not related to alpha from Likelihood
        alpha_prop = alpha * prop[idx_g]
        return -(
            n_cells * loggamma(alpha)
            - np.sum(loggamma(tb + alpha))
            + np.sum(loggamma(y + alpha_prop))
            - np.sum(loggamma(alpha_prop))
        )

    estimator = so.minimize_scalar(neg_logL, bounds=bounds, method="bounded")
    return estimator.x
