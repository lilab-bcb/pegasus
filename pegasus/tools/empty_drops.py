import itertools
import numba

import numpy as np
import pandas as pd
import scipy.optimize as so

from pegasusio import MultimodalData, UnimodalData
from pegasus.tools import eff_n_jobs
from pegasus.tools import SimpleGoodTuring
from pegasus.cylib.fast_utils import test_empty_drops

from anndata import AnnData
from scipy.special import loggamma
from scipy.sparse import csr_matrix, issparse
from scipy.interpolate import UnivariateSpline
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
    alpha: Optional[float] = None,
    n_iters: int = 10000,
    random_state: int = 0,
    exclude_from: int = 50,
    significance_level: float = 0.05,
    n_jobs: int = -1,
    ambient_proportion: Optional[np.array] = None,  # Testing
    knee: Optional[float] = None,  # Testing
) -> None:
    if "counts" in data._unidata.matrices:
        mat_key = "counts"
    elif "raw.X" in data._unidata.matrices:
        mat_key = "raw.X"
    X = data.get_matrix(mat_key) if mat_key else data.X

    if "n_counts" not in data.obs:
        logger.info(
            "Calculate n_counts for EmptyDrops since 'n_counts' not in data.obs."
        )
        df_droplets = pd.DataFrame(index=data.obs_names)
        df_droplets["n_counts"] = X.sum(axis=1).A1
    else:
        df_droplets = data.obs[["n_counts"]].copy()

    assert (
        np.sum(X.sum(axis=0) == 0) == 0
    ), "Some genes have zero counts in data! Please run pegasus.identify_robust_genes() first!"

    idx_low = df_droplets.index.get_indexer(
        df_droplets.loc[df_droplets["n_counts"] <= thresh_low].index.values
    )
    G = X[idx_low, :]

    if ambient_proportion is None:
        ambient_profile = G.sum(axis=0).A1
        sgt = SimpleGoodTuring(ambient_profile)
        ambient_proportion = sgt.get_proportions()

    if alpha is None:
        alpha = _estimate_alpha(G, ambient_proportion)
        logger.info(f"Calculating alpha is finished. Estimate alpha = {alpha}.")

    cells_test = df_droplets.loc[df_droplets["n_counts"] > thresh_low].index.values
    idx_test = df_droplets.index.get_indexer(cells_test)
    logger.info(f"Run test through {idx_test.size} cells with n_counts > {thresh_low}.")

    T = X[idx_test, :]
    tb = T.sum(axis=1).A1 if issparse(T) else T.sum(axis=1)
    Lb_data = _logL_data_dep(T, ambient_proportion, alpha)
    Lb_alpha = _logL_data_indep(tb, alpha)

    n_jobs = eff_n_jobs(n_jobs)
    n_below = _test_empty_drops(alpha, ambient_proportion, tb, Lb_data, n_iters, random_state, n_jobs)
    pval = (n_below + 1) / (n_iters + 1)
    logger.info("Calculation on p-values is finished.")

    df_rank_stats, knee_est, inflection = _rank_barcode(X, thresh_low, exclude_from)
    if knee is None:
        knee = knee_est

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

    return (df_rank_stats, knee, inflection)


def _logL_data_dep(Y, prop, alpha):
    if issparse(Y):
        idx_b, idx_g = Y.nonzero()
        y = Y.data
    else:
        idx_b, idx_g = np.nonzero(Y)
        y = Y[idx_b, idx_g]
    alpha_prop = alpha * prop[idx_g]
    L_a1 = loggamma(y + alpha_prop) - loggamma(y + 1) - loggamma(alpha_prop)
    L_mat = csr_matrix((L_a1, (idx_b, idx_g)), shape=Y.shape)
    return L_mat.sum(axis=1).A1


def _logL_data_indep(t_b, alpha):
    return loggamma(t_b + 1) + loggamma(alpha) - loggamma(t_b + alpha)


def _rank_barcode(X, thresh_low, exclude_from):
    total_counts = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
    rank_values, rank_sizes = np.unique(total_counts, return_counts=True)
    rank_values = rank_values[::-1]
    rank_sizes = rank_sizes[::-1]
    tie_rank = np.cumsum(rank_sizes) - (rank_sizes - 1) / 2  # Get mid-rank of each tie

    idx_keep = np.where(rank_values > thresh_low)[0]
    y = np.log1p(rank_values[idx_keep])
    x = np.log(tie_rank[idx_keep])
    left_edge, right_edge = _find_curve_bounds(x, y, exclude_from)

    inflection = np.expm1(y[right_edge])

    fitted_values = np.log1p(rank_values)
    idx_focus = np.arange(left_edge, right_edge+1)
    if idx_focus.size >= 4:
        curx = x[idx_focus]
        cury = y[idx_focus]
        xbounds = curx[[0, -1]]
        ybounds = cury[[0, -1]]
        slope = (ybounds[1] - ybounds[0]) / (xbounds[1] - xbounds[0])
        intercept = ybounds[0] - xbounds[0] * slope
        above = np.where(cury >= curx * slope + intercept)[0]
        dist = np.abs(slope * curx[above] - cury[above] + intercept) / np.sqrt(slope**2 + 1)
        knee = np.expm1(cury[above[np.argmax(dist)]])
    else:
        knee = np.expm1(y[idx_focus[0]])

    def repeat(vals, lens):
        return [element for element, count in zip(vals, lens) for _ in range(count)]

    df_stats = pd.DataFrame({
        "Barcodes": repeat(np.log(tie_rank), rank_sizes),
        "n_counts": repeat(rank_values, rank_sizes),
        "UMI counts": repeat(fitted_values, rank_sizes),
    })
    return (df_stats, knee, inflection)


def _find_curve_bounds(x, y, exclude_from):
    d1n = np.diff(y) / np.diff(x)

    skip = np.min([d1n.size - 1, np.sum(x<=np.log(exclude_from))])
    d1n = d1n[-(d1n.size - skip):]

    right_edge = np.argmin(d1n)
    left_edge = np.argmax(d1n[:(right_edge+1)])

    return (left_edge + skip, right_edge + skip)


def _test_empty_drops(alpha, prop, tb, P_data, n_iters, random_state, n_jobs, temp_folder=None):
    idx_sorted = np.lexsort((P_data, tb))
    P_sorted = P_data[idx_sorted]
    tb_unique, tb_cnt = np.unique(tb[idx_sorted], return_counts=True)

    alpha_prop = alpha * prop
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
