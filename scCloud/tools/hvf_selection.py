import time
import numpy as np
import numpy.ma as ma
import pandas as pd

from scipy.sparse import issparse
from collections import defaultdict
from joblib import Parallel, delayed
import skmisc.loess as sl
from typing import List



def calc_ba_mean_and_var(data: 'AnnData') -> None:
    """ Calculate batch adjusted mean and variance
    """

    assert ('Channels' in data.uns) and ('Groups' in data.uns)
    assert issparse(data.X)

    channels = data.uns['Channels']
    group_dict = defaultdict(list)

    means = np.zeros((data.shape[1], channels.size))
    partial_sum = np.zeros((data.shape[1], channels.size))
    ncells = np.zeros(channels.size)

    for i, channel in enumerate(channels):
        idx = np.isin(data.obs['Channel'], channel)
        mat = data.X[idx].astype(np.float64)

        ncells[i] = mat.shape[0]
        if ncells[i] > 0:
            if ncells[i] == 1:
                means[:, i] = mat.toarray()[0]
            else:
                means[:, i] = mat.mean(axis=0).A1
                m2 = mat.power(2).sum(axis=0).A1
                partial_sum[:, i] = m2 - ncells[i] * (means[:, i] ** 2)

            group = data.obs['Group'][idx.nonzero()[0][0]]
            group_dict[group].append(i)

    partial_sum[partial_sum < 1e-6] = 0.0

    overall_means = np.dot(means, ncells) / data.shape[0]
    batch_adjusted_vars = np.zeros(data.shape[1])

    groups = data.uns['Groups']
    for group in groups:
        gchannels = group_dict[group]
        ncell = ncells[gchannels].sum()
        gm = np.dot(means[:, gchannels], ncells[gchannels]) / ncell

        if groups.size > 1:
            batch_adjusted_vars += ncell * ((gm - overall_means) ** 2)

    data.var['mean'] = overall_means
    data.var['var'] = (batch_adjusted_vars + partial_sum.sum(axis=1)) / (
        data.shape[0] - 1.0
    )


def select_hvf_scCloud(
    data: 'AnnData', consider_batch: bool, n_top: int = 2000, span: float = 0.02, benchmark_time: bool = False
) -> None:
    """ Select highly variable features using the scCloud method
    """
    if 'robust' not in data.var:
        raise ValueError("Please run `qc_metrics` to identify robust genes")
    robust_idx = data.var['robust'].values
    hvf_index = np.zeros(robust_idx.sum(), dtype=bool)

    if consider_batch:
        if benchmark_time or ('mean' not in data.var) or ('var' not in data.var):
            calc_ba_mean_and_var(data)
        mean = data.var.loc[robust_idx, 'mean']
        var = data.var.loc[robust_idx, 'var']
    else:
        X = data.X[:, robust_idx]
        mean = X.mean(axis=0).A1
        m2 = X.power(2).sum(axis=0).A1
        var = (m2 - X.shape[0] * (mean ** 2)) / (X.shape[0] - 1)

        data.var['mean'] = 0.0
        data.var['var'] = 0.0
        data.var.loc[robust_idx, 'mean'] = mean
        data.var.loc[robust_idx, 'var'] = var

    lobj = sl.loess(mean, var, span=span, degree=2)
    lobj.fit()

    rank1 = np.zeros(hvf_index.size, dtype=int)
    rank2 = np.zeros(hvf_index.size, dtype=int)

    delta = var - lobj.outputs.fitted_values
    fc = var / lobj.outputs.fitted_values

    rank1[np.argsort(delta)[::-1]] = range(hvf_index.size)
    rank2[np.argsort(fc)[::-1]] = range(hvf_index.size)
    hvf_rank = rank1 + rank2

    hvf_index[np.argsort(hvf_rank)[:n_top]] = True

    data.var['hvf_loess'] = 0.0
    data.var.loc[robust_idx, 'hvf_loess'] = lobj.outputs.fitted_values

    data.var['hvf_rank'] = -1
    data.var.loc[robust_idx, 'hvf_rank'] = hvf_rank
    data.var['highly_variable_features'] = False
    data.var.loc[robust_idx, 'highly_variable_features'] = hvf_index


def select_hvf_seurat_single(
    X: 'csr_matrix', n_top: int, min_disp: float, max_disp: float, min_mean: float, max_mean: float) -> List[int] :
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

    df = pd.DataFrame({'log_dispersion': dispersion, 'bin': pd.cut(mean, bins=20)})
    log_disp_groups = df.groupby('bin')['log_dispersion']
    log_disp_mean = log_disp_groups.mean()
    log_disp_std = log_disp_groups.std(ddof=1)
    log_disp_zscore = (
        df['log_dispersion'].values - log_disp_mean.loc[df['bin']].values
    ) / log_disp_std.loc[df['bin']].values
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


def select_hvf_seurat_multi(X: 'csr_matrix', channels: List[str], cell2channel: List[str], n_top: int, n_jobs: int) -> List[int]:
    Xs = []
    for channel in channels:
        Xs.append(X[np.isin(cell2channel, channel)])

    from joblib import effective_n_jobs
    n_jobs = effective_n_jobs(n_jobs)

    res_arr = np.array(
        Parallel(n_jobs=n_jobs)(
            delayed(select_hvf_seurat_single)(Xs[i], n_top)
            for i in range(channels.size)
        )
    )
    selected = res_arr >= 0
    shared = selected.sum(axis=0)
    cands = (shared > 0).nonzero()[0]
    median_rank = ma.median(ma.masked_array(res_arr, mask=~selected), axis=0).data
    cands = sorted(cands, key=lambda x: median_rank[x])
    cands = sorted(cands, key=lambda x: shared[x], reverse=True)

    hvf_rank = np.full(X.shape[1], -1, dtype=int)
    hvf_rank[cands[:n_top]] = range(n_top)

    return hvf_rank


def select_hvf_seurat(
    data: 'AnnData',
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

    robust_idx = data.var['robust'].values
    X = data.X[:, robust_idx]

    hvf_rank = (
        select_hvf_seurat_multi(
            X, data.uns['Channels'], data.obs['Channel'], n_top, n_jobs=n_jobs
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

    data.var['hvf_rank'] = -1
    data.var.loc[robust_idx, 'hvf_rank'] = hvf_rank
    data.var['highly_variable_features'] = False
    data.var.loc[robust_idx, 'highly_variable_features'] = hvf_index


def highly_variable_features(
    data: 'AnnData',
    consider_batch: bool,
    flavor: str = "scCloud",
    n_top: int = 2000,
    span: float = 0.02,
    min_disp: float = 0.5,
    max_disp: float = np.inf,
    min_mean: float = 0.0125,
    max_mean: float = 7,
    n_jobs: int = -1,
    benchmark_time: bool =False,
) -> None:
    """
    TODO: Documentation.
    n_top == None refers to no n_top genes
    n_jobs: -1 refers to use all available threads
    """

    start = time.time()

    if consider_batch and 'Channels' not in data.uns:
        print(
            "Warning: Batch correction parameters are not calculated. Switch to not considering batch for variable gene selection."
        )
        consider_batch = False

    if flavor == 'scCloud':
        select_hvg_scCloud(
            data,
            consider_batch,
            n_top=n_top,
            span=span,
            benchmark_time=benchmark_time,
        )
    else:
        assert flavor == 'Seurat'
        select_hvg_seurat(
            data,
            consider_batch,
            n_top=n_top,
            min_disp=min_disp,
            max_disp=max_disp,
            min_mean=min_mean,
            max_mean=max_mean,
            n_jobs=n_jobs,
        )

    end = time.time()
    print(
        "{tot} highly variable features have been selected. Time spent = {time:.2f}s.".format(tot = data.var['highly_variable_features'].sum(), time = end - start)
    )


def collect_highly_variable_gene_matrix(data):
    data_c = data[:, data.var['highly_variable_genes']].copy()
    data_c.X = data_c.X.toarray()

    return data_c
