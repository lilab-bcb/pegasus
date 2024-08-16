import itertools

import numpy as np
import pandas as pd
import scipy.optimize as so

from pegasusio import MultimodalData, UnimodalData
from pegasus.tools import eff_n_jobs
from anndata import AnnData
from scipy.special import loggamma
from scipy.sparse import coo_matrix
from scipy.interpolate import UnivariateSpline
from statsmodels.stats.multitest import fdrcorrection as fdr
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
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
    chunk_size: int = 5000,
    n_jobs: int = -1,
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
    ambient_profile = G.sum(axis=0).A1
    ambient_proportion = _simple_good_turing(ambient_profile)

    if not alpha:
        alpha = _estimate_alpha(G, ambient_proportion)
        logger.info(f"Calculating alpha is finished. Estimate alpha = {alpha}.")

    t_b = X.sum(axis=1).A1
    L_b = _logL(alpha, ambient_proportion, t_b, X)
    n_cells = t_b.size

    n_jobs = eff_n_jobs(n_jobs)
    n_below = _test_ambient(alpha, ambient_proportion, t_b, L_b, n_iters, random_state, chunk_size, n_jobs)
    logger.info("Calculation on p-values is finished.")
    pval = (n_below + 1) / (n_iters + 1)
    _, qval = fdr(pval)

    data.obs["EmptyDrops.pval"] = pval
    data.obs["EmptyDrops.qval"] = qval

    df_rank_stats, knee, inflection = _rank_barcode(t_b, thresh_low, exclude_from)
    logger.info("Calculate the knee point.")

    data.obs["EmptyDrops.quality"] = "low"
    data.obs.loc[data.obs["EmptyDrops.qval"]<significance_level, "EmptyDrops.quality"] = "high"
    data.obs.loc[data.obs["n_counts"]>inflection, "EmptyDrops.quality"] = "high"


def _rank_barcode(total_counts, thresh_low, exclude_from):
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
        spline = UnivariateSpline(x[idx_focus], y[idx_focus], k=3, s=5)
        fitted_values[idx_focus] = spline(x[idx_focus])
        d1 = spline.derivative(n=1)(x[idx_focus])
        d2 = spline.derivative(n=2)(x[idx_focus])
        curvature = d2 / (1 + d1**2)**1.5
        knee = np.expm1(y[idx_focus][np.argmin(curvature)])
    else:
        knee = np.expm1(y[idx_focus[0]])

    def repeat(vals, lens):
        return [element for element, count in zip(vals, lens) for _ in range(count)]

    df_stats = pd.DataFrame({
        "Barcodes": repeat(tie_rank, rank_sizes),
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


def _test_ambient(alpha, prop, t_b, L_b, n_iters, random_state, chunk_size, n_jobs):
    n_cells = t_b.size
    n_genes = prop.size
    n_below = np.zeros(n_cells, dtype=int)
    np.random.seed(random_state)

    p_arr = np.random.dirichlet(alpha * prop, size=n_iters)
    print(p_arr.shape)
    chunks = process_map(
        _test_ambient_by_chunk,
        itertools.repeat(alpha),
        itertools.repeat(prop),
        itertools.repeat(p_arr),
        [t_b[i:(i+chunk_size)] for i in range(0, n_cells, chunk_size)],
        [L_b[i:(i+chunk_size)] for i in range(0, n_cells, chunk_size)],
        tqdm_class=tqdm,
        max_workers=n_jobs,
    )
    return np.hstack(chunks)


def _test_ambient_by_chunk(alpha, prop, p_arr, t_b, L_b):
    n_below = np.zeros(t_b.size, dtype=int)
    for p in p_arr:
        samples = coo_matrix([np.random.multinomial(t, p) for t in t_b])
        L_test = _logL(alpha, prop, t_b, samples)
        n_below += (L_test <= L_b)
    return n_below


def _simple_good_turing(counts):
    from nltk.probability import FreqDist, SimpleGoodTuringProbDist
    from collections import Counter

    n_genes = counts.size
    idx_obs = np.where(counts > 0)[0]
    n0 = n_genes - idx_obs.size
    freq_dist = FreqDist(Counter(dict(zip(np.arange(n_genes), counts))))
    sgt = SimpleGoodTuringProbDist(freq_dist)

    ambient_proportion = (
        sgt.prob(None) / n0 * np.ones(n_genes) if n0 != 0 else np.zeros(n_genes)
    )
    for i in idx_obs:
        ambient_proportion[i] = sgt.prob(i)

    return ambient_proportion


def _estimate_alpha(G, prop):
    total_counts = G.sum(axis=0).A1
    n_cells = total_counts.size
    indices = G.nonzero()[1]
    x = G.data

    def neg_logL(alpha):
        return -(
            np.sum(loggamma(total_counts + 1))
            + n_cells * loggamma(alpha)
            - np.sum(loggamma(total_counts + alpha))
            + np.sum(loggamma(x + alpha * prop[indices]))
            - np.sum(loggamma(x + 1))
            - np.sum(n_cells * loggamma(alpha * prop))
        )

    estimator = so.minimize(neg_logL, np.array([1]), bounds=so.Bounds(0.01, 10000))
    return estimator.x[0]


def _logL(alpha, prop, total_counts, Y):
    n_cells = total_counts.size
    n_genes = prop.size
    idx_cells, idx_genes = Y.nonzero()
    x = Y.data
    Y_ag = coo_matrix((loggamma(x + alpha * prop[idx_genes]), (idx_cells, idx_genes)), shape=(n_cells, n_genes))
    Y_p1 = coo_matrix((loggamma(x + 1), (idx_cells, idx_genes)), shape=(n_cells, n_genes))

    return (
        loggamma(total_counts + 1)
        + loggamma(alpha)
        - loggamma(total_counts + alpha)
        + np.asarray(Y_ag.sum(axis=1)).squeeze()
        - np.asarray(Y_p1.sum(axis=1)).squeeze()
        - np.sum(loggamma(alpha * prop))
    )
