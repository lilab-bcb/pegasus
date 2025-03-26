import numpy as np
import pandas as pd
import scipy.optimize as so

from pegasusio import MultimodalData, UnimodalData
from pegasus.tools import sgt_estimate
from pegasus.plotting import plot_barcode_rank

from numba import njit
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
def empty_drops_numba(
    data: Union[MultimodalData, UnimodalData, AnnData],
    mat_key: Optional[str] = None,
    thresh_low: int = 100,
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
    """Estimation on cell barcodes to decide if ambient or real cells by EmptyDrops algorithm. [Lun19]_

    Parameters
    -----------
    data: ``pegasusio.MultimodalData`` or ``pegasusio.UnimodalData`` or ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    mat_key: ``str``, optional, default: ``None``
        The count matrix key to use. By default use ``data.X``.

    thresh_low: ``int``, optional, default ``100``
        The threshold to determine the ambient profile.
        Cells with UMI counts <= ``thresh_low`` are assigned to ambient initially.

    distribution: ``str``, optinal, default: ``multi-dirichlet``
        The distribution of UMI counts across genes for ambient profile, i.e. Null Hypothesis.
        By default, assume Multinomial Dirichlet. Alternatively, set ``"multinomial"`` to assume Multinomial distribution (i.e. without ``alpha``).

    alpha: `float``, optional, default: ``None``
        The alpha term used for defining the distribution. Only take effects when ``distribution="multi-dirichlet"``.
        By default, estimate from data. Otherwise, users can specify custom ``alpha`` which needs to be positive.

    n_iters: ``int``, optional, default: ``10000``
        Number of iterations of Monte-Carlo simulation for calculating p-values.

    significance_level: ``float``, optional, default: ``0.01``
        Significance level applied to FDR q-values to decide non-ambient cells.

    random_state: ``int``, optional, default: ``0``
        Random state for reproducibility.

    knee_method: ``str``, optional, default: ``"linear"``
        Regression method used to find the knee point.
        By default, use linear regression. Alternatively, users can set to ``"spline"` to use the smooth spline method which is described in paper.

    exclude_from: ``int``, optional, default: ``50``
        Parameter used for finding the knee point. Specify the number of highest ranked cell barcodes to exclude from fitting so as to avoid problems with discreteness.

    show_summary_plot: ``bool``, optional, default: ``False``
        Whether plotting barcode rank plot or not. By default ``False``.

    n_jobs: ``int``, optional, default: ``1``
        Number of jobs used in Monte-Carlo simulation.
        Default is 1 job. Specify ``-1`` to use all available vCPUs.

    Returns
    -------
    Knee point and inflection point calculated.

    Update ``data.obs``:
        * ``data.obs["ambient"]``: Boolean attribute. ``True`` if the cell is in ambient; ``False`` if it's the real cell.
        * ``data.obs["EmptyDrops.logL"]``: The log-likelihood of each cell.
        * ``data.obs["EmptyDrops.pval"]``: p-values of cells from the EmptyDrops hypothesis test.
        * ``data.obs["EmptyDrops.qval"]``: FDR q-values of cells from the EmptyDrops  hypothesis test.
    """
    if mat_key is None:
        if "counts" in data._unidata.matrices:
            mat_key = "counts"
        elif "raw.X" in data._unidata.matrices:
            mat_key = "raw.X"
    X = data.get_matrix(mat_key) if mat_key else data.X

    assert (
        np.sum(X.sum(axis=0) == 0) == 0
    ), "Some genes have zero counts in data! Please run pegasus.identify_robust_genes() first!"

    if "n_counts" not in data.obs:
        logger.info(
            "Calculate n_counts for EmptyDrops since 'n_counts' not in data.obs."
        )
        n_counts = X.sum(axis=1).A1 if issparse(X) else X.sum(axis=1)
    else:
        n_counts = data.obs["n_counts"].values

    idx_low = np.where(n_counts <= thresh_low)[0]
    G = X[idx_low, :]

    ambient_profile = G.sum(axis=0).A1
    ambient_proportion = sgt_estimate(ambient_profile)

    assert distribution in ["multi-dirichlet", "multinomial"], "distribution must be chosen from ['multi-dirichlet', 'multinomial']!"
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

    n_below = _test_empty_drops_numba(use_alpha, alpha, ambient_proportion, tb, Lb_data, n_iters, random_state, n_jobs)
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

    return (knee, inflection)


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


def _logL_data_dep(Y, prop, alpha, use_alpha):
    """Calculate the count matrix dependent part of logL, i.e. terms within the summation notation.
    """
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
    """Calculate the count matrix independent part of logL, i.e. terms outside the summation notation.
    """
    if use_alpha:
        return loggamma(tb + 1) + loggamma(alpha) - loggamma(tb + alpha)
    else:
        return loggamma(tb + 1)


def _test_empty_drops_numba(use_alpha, alpha, prop, tb, P_data, n_iters, random_state, n_jobs, temp_folder=None):
    # Sort cells by UMI counts in ascending order, then for ties sort by logL in ascending order
    idx_sorted = np.lexsort((P_data, tb))
    P_sorted = P_data[idx_sorted]
    tb_unique, tb_cnt = np.unique(tb[idx_sorted], return_counts=True)

    # Debug
    #gs_dbg = np.loadtxt("/home/ubuntu/empty_drops_debug/GS_R.txt").astype(int)

    # Parameters passed to Monte Carlo simulation
    alpha_prop = alpha * prop if use_alpha else prop
    n_cells = tb.size

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
            # Monte-Carlo simulation
            delayed(_test_by_chunk)(
                use_alpha,
                alpha_prop,
                tb_unique,
                tb_cnt,
                P_sorted,
                n_cells,
                random_state,   #TODO: Only 1 random state provided here, as my test cases always have n_jobs=1; modify if final decision is to keep this parallel feature
                intervals[i][0],
                intervals[i][1],
                #gs_dbg,
            )
            for i in range(n_jobs)
        )
    logger.info("Significance test is finished.")
    n_below_sorted = np.vstack(result).sum(axis=0)

    # Get n_below (ordered by original cell order) from n_below_sorted (ordered by P_sorted)
    idx_inv = np.empty_like(idx_sorted)
    idx_inv[idx_sorted] = np.arange(idx_sorted.size)
    n_below = n_below_sorted[idx_inv]

    return n_below


@njit
def _test_by_chunk(
    use_alpha: bool,
    alpha_prop: np.array,
    tb_unique: np.array,
    tb_cnt: np.array,
    L_b: np.array,
    n_cells: int,
    random_state: int,
    start_pos: int,
    end_pos: int,
    #gs_dbg: np.array,
) -> np.array:
    tb_max = tb_unique[-1]
    n_genes = alpha_prop.size

    below_cnt = np.zeros(n_cells, dtype=np.int64)

    np.random.seed(random_state)
    for i in np.arange(start_pos, end_pos):
        prob_arr = dirichlet_sampling(alpha_prop) if use_alpha else alpha_prop
        gs_arr = weighted_sampling(prob_arr, tb_max)
        #gs_arr = gs_dbg[i, :]

        z = np.zeros(n_genes, dtype=np.int64)
        idx_tb = 0
        idx_b = 0
        idx_g = 0
        cur_L = 0.0
        t = 0

        while idx_tb < tb_unique.size:
            while t < tb_unique[idx_tb]:
                k = gs_arr[idx_g]
                z[k] += 1
                if use_alpha:
                    cur_L += (np.log(z[k] + alpha_prop[k] - 1) - np.log(z[k]))
                else:
                    cur_L += (np.log(alpha_prop[k]) - np.log(z[k]))
                t += 1
                idx_g += 1

            # First idx_below satisfying L_b[idx_below - 1] < cur_L <= L_b[idx_below]
            idx_below = np.searchsorted(L_b[idx_b:(idx_b+tb_cnt[idx_tb])], cur_L, side="left")
            if idx_below < tb_cnt[idx_tb]:
                below_cnt[idx_b+idx_below] += 1

            idx_b += tb_cnt[idx_tb]
            idx_tb += 1

        assert idx_g == tb_max, f"Only {idx_g} counts out of {tb_max} are used!"

    # Use the same way as R package, which returns n_chunk vectors of shape (n_cells,).
    # While in Cell Ranger implementation, Monte-Carlo step returns 2d vector of shape (tb_unique.size, n_iters), which takes more space.
    idx_below = 0
    for idx_tb in np.arange(tb_unique.size):
        for i in np.arange(tb_cnt[idx_tb] - 1):
            prev = below_cnt[idx_below]
            idx_below += 1
            below_cnt[idx_below] += prev
        idx_below += 1

    return below_cnt


@njit
def weighted_sampling(weights, n_samples):
    prob_cum = np.cumsum(weights)
    samples = np.empty(n_samples, dtype=np.int64)
    for i in np.arange(n_samples):
        rand_val = np.random.random()
        # Last j satisfying prob_cum[j-1] <= rand_val < prob_cum[j]
        j = np.searchsorted(prob_cum, rand_val, side="right")
        samples[i] = j
    return samples


@njit
def dirichlet_sampling(alpha_prop):
    res = np.full(alpha_prop.size, np.nan)
    s = 0.0
    for i in np.arange(alpha_prop.size):
        res[i] = np.random.gamma(shape=alpha_prop[i], scale=1.0)
        s += res[i]
    return res / s


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
    else:
        knee = int(np.ceil(np.exp(y[idx_focus[0]])))

    df_rank = pd.DataFrame({
        "rank": tie_rank,
        "size": rank_sizes,
        "n_counts_obs": rank_values,
        "n_counts_fitted": fitted,
    })

    return (df_rank, knee, inflection)


def _find_knee_point(x, y, fitted, idx_focus, method="linear"):
    x_obs = x[idx_focus]
    y_obs = y[idx_focus]

    if method == "spline":  # Described in paper
        from scipy.interpolate import UnivariateSpline

        spline = UnivariateSpline(x_obs, y_obs, k=3, s=5)    # Try to approximate smooth.spline() in R, which only uses df parameter.
        fitted[idx_focus] = spline(x_obs)
        d1 = spline.derivative(n=1)(x_obs)
        d2 = spline.derivative(n=2)(x_obs)
        curvature = d2 / (1 + d1**2)**1.5
        knee = int(np.ceil(np.exp(y_obs[np.argmin(curvature)])))
    else:    # Used in R package's latest version.
        slope = (y_obs[-1] - y_obs[0]) / (x_obs[-1] - x_obs[0])
        intercept = y_obs[0] - x_obs[0] * slope
        y_fitted = x_obs * slope + intercept
        fitted[idx_focus] = y_fitted
        above = np.where(y_obs >= y_fitted)[0]
        distance = (y_obs[above] - y_fitted[above]) / np.sqrt(slope**2 + 1)
        knee = int(np.ceil(np.exp(y_obs[above[np.argmax(distance)]])))

    return knee


@njit
def _find_curve_bounds(x, y, exclude_from):
    d1n = np.diff(y) / np.diff(x)

    skip = np.min(np.array([d1n.size - 1, np.sum(x<=np.log(exclude_from))]))
    d1n = d1n[-(d1n.size - skip):]

    right_edge = np.argmin(d1n)   # point with the least steep
    left_edge = np.argmax(d1n[:(right_edge+1)])   # point with the highest steep to the left of right edge

    return (left_edge + skip, right_edge + skip)
