import itertools
import numba

import numpy as np
import pandas as pd
import scipy.optimize as so

from pegasusio import MultimodalData, UnimodalData
from pegasus.tools import eff_n_jobs
from pegasus.tools import SimpleGoodTuring
from anndata import AnnData
from scipy.special import loggamma, factorial
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
    sgt = SimpleGoodTuring(ambient_profile)
    ambient_proportion = sgt.get_proportions()

    if not alpha:
        alpha = _estimate_alpha(G, ambient_proportion)
        logger.info(f"Calculating alpha is finished. Estimate alpha = {alpha}.")

    cells_test = df_droplets.loc[df_droplets["n_counts"] > thresh_low].index.values
    idx_test = df_droplets.index.get_indexer(cells_test)
    logger.info(f"Run test through {idx_test.size} cells passing thresh_low.")

    T = X[idx_test, :]
    t_b = T.sum(axis=1).A1
    L_b = _logL(alpha, ambient_proportion, t_b, T)

    n_jobs = eff_n_jobs(n_jobs)
    n_below = _test_empty_drops(alpha, ambient_proportion, t_b, L_b, n_iters, random_state, n_jobs)
    pval = (n_below + 1) / (n_iters + 1)
    logger.info("Calculation on p-values is finished.")

    df_rank_stats, knee, inflection = _rank_barcode(t_b, thresh_low, exclude_from)
    idx_always = np.where(t_b >= knee)[0]
    pval[idx_always] = 0.0
    logger.info("Adjust p-values by the estimated knee point.")

    _, qval = fdr(pval)
    logger.info

    data.obs["EmptyDrops.pval"] = np.nan
    data.obs.loc[cells_test, "EmptyDrops.pval"] = pval
    data.obs["EmptyDrops.qval"] = np.nan
    data.obs.loc[cells_test, "EmptyDrops.qval"] = qval

    return (df_rank_stats, knee, inflection)


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


def _test_empty_drops(alpha, prop, t_b, L_b, n_iters, random_state, n_jobs):
    tb_dict = {}
    for idx, t in enumerate(t_b):
        if t in tb_dict:
            tb_dict[t].append(idx)
        else:
            tb_dict[t] = [idx]
    tb_unique = list(tb_dict.keys())
    tb_unique.sort()

    np.random.seed(random_state)
    p_array = np.random.dirichlet(alpha * prop, size=n_iters)
    chunk_size = n_iters // n_jobs
    chunks = process_map(
        _test_empty_drops_by_chunk,
        itertools.repeat(alpha),
        itertools.repeat(prop),
        itertools.repeat(tb_unique),
        itertools.repeat(tb_dict),
        itertools.repeat(t_b),
        itertools.repeat(L_b),
        itertools.repeat(random_state),
        [p_array[i:(i+chunk_size)] for i in range(0, n_iters, chunk_size)],
        tqdm_class=tqdm,
        max_workers=n_jobs,
    )
    return np.vstack(chunks).sum(axis=0)


def _test_empty_drops_by_chunk(alpha, prop, tb_unique, tb_dict, t_b, L_b, random_state, p_arr):
    tb_min = tb_unique[0]
    idx_all_genes = np.arange(prop.size)
    n_below = np.zeros(t_b.size)
    L_test = np.zeros(t_b.size)

    np.random.seed(random_state)
    for p in p_arr:
        for t in tb_unique:
            if t == tb_min:
                z = np.random.multinomial(t, p)
                L = _logL(alpha, prop, np.array(t), coo_matrix(z))
            else:
                idx_extra_genes = np.random.choice(idx_all_genes, size=t-t1, replace=True, p=p)
                for i in idx_extra_genes:
                    z[i] += 1
                L = _calc_logL_by_increment(L1, alpha, prop, t, t - t1, z, idx_extra_genes)

            L_test[tb_dict[t]] = L
            t1 = t
            L1 = L
        n_below += (L_test <= L_b)

    return n_below


def _calc_logL_by_increment(L1, alpha, prop, t, n_extra_counts, z, idx_extra_genes):
    return (
        L1
        + n_extra_counts * np.log(t / (t + alpha - 1))
        + np.sum(np.log((z[idx_extra_genes] + alpha * prop[idx_extra_genes] - 1) / z[idx_extra_genes]))
    )


#def _test_empty_drops(alpha, prop, t_b, L_b, n_iters, random_state, chunk_size, n_jobs):
#    n_cells = t_b.size
#    n_genes = prop.size
#    n_below = np.zeros(n_cells, dtype=int)
#    np.random.seed(random_state)
#
#    #p = np.random.dirichlet(alpha * prop, size=n_iters)
#    n_below = np.zeros(n_cells, dtype=int)
#    for k in range(n_iters):
#        p = np.random.dirichlet(alpha * prop)
#        chunks = process_map(
#            _test_ambient_by_chunk,
#            itertools.repeat(alpha),
#            itertools.repeat(prop),
#            itertools.repeat(p),
#            [t_b[i:(i+chunk_size)] for i in range(0, n_cells, chunk_size)],
#            [L_b[i:(i+chunk_size)] for i in range(0, n_cells, chunk_size)],
#            tqdm_class=tqdm,
#            max_workers=n_jobs,
#        )
#        n_below += np.hstack(chunks)
#        if (k + 1) % 1000 == 0:
#            print(f"{k+1} iterations finished.")
#    return n_below


#def _test_ambient_by_chunk(alpha, prop, p, t_b, L_b):
#    n_below = np.zeros(t_b.size, dtype=int)
#    samples = coo_matrix([np.random.multinomial(t, p) for t in t_b])
#    L_test = _logL(alpha, prop, t_b, samples)
#    return (L_test <= L_b)


def _estimate_alpha(G, prop):
    total_counts = G.sum(axis=1).A1
    n_cells = total_counts.size
    indices = G.nonzero()[1]
    x = G.data

    def neg_logL(alpha):
        # TODO: Maybe remove terms not related to alpha
        return -(
            np.sum(loggamma(total_counts + 1))
            + n_cells * loggamma(alpha)
            - np.sum(loggamma(total_counts + alpha))
            + np.sum(loggamma(x + alpha * prop[indices]))
            - np.sum(loggamma(x + 1))
            - np.sum(n_cells * loggamma(alpha * prop[indices]))
        )

    estimator = so.minimize(neg_logL, np.array([1]), bounds=so.Bounds(0.01, 10000))
    return estimator.x[0]


def _logL(alpha, prop, total_counts, Y):
    n_cells = total_counts.size
    n_genes = prop.size
    idx_cells, idx_genes = Y.nonzero()
    x = Y.data
    Y_ag = coo_matrix((loggamma(x + alpha * prop[idx_genes]), (idx_cells, idx_genes)), shape=(n_cells, n_genes))
    Y_fact = coo_matrix((loggamma(x + 1), (idx_cells, idx_genes)), shape=(n_cells, n_genes))
    A_g = coo_matrix((loggamma(alpha * prop[idx_genes]), (idx_cells, idx_genes)), shape=(n_cells, n_genes))

    return (
        loggamma(total_counts + 1)
        + loggamma(alpha)
        - loggamma(total_counts + alpha)
        + np.asarray(Y_ag.sum(axis=1)).squeeze()
        - np.asarray(Y_fact.sum(axis=1)).squeeze()
        - np.asarray(A_g.sum(axis=1)).squeeze()
    )
