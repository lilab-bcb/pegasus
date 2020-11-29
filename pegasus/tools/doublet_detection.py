import time
import numpy as np
import pandas as pd

from typing import List, Optional, Union, Tuple

import logging
logger = logging.getLogger(__name__)

from pegasusio import UnimodalData, MultimodalData
from pegasusio import timer



def _f1(f, x, h): # calculated using five-point stencil
    if x - 2 < 0 or x + 2 >= f.size:
        return np.nan
    return (-f[x + 2] + 8 * f[x + 1] - 8 * f[x - 1] + f[x - 2]) / 12 / h

def _f2(f, x, h): # calculated using five-point stencil
    if x - 2 < 0 or x + 2 >= f.size:
        return np.nan
    return (-f[x + 2] + 16 * f[x + 1] - 30 * f[x] + 16 * f[x - 1] - f[x - 2]) / 12 / h / h

def _curvature(f, x, h): # calculated curvature
    return _f2(f, x, h) / (1.0 + _f1(f, x, h) ** 2) ** 1.5

def _calc_vec_f(func, size, f, h): # convenient function to vetorize the above functions
    res = np.zeros(size)
    for i in range(size):
        res[i] = func(f, i, h)
    return res

def _find_local_maxima(y: List[float], frac: float = 1.0 / 3.0, merge_peak_frac: float = 0.06) -> Tuple[List[int], List[int], List[int]]:
    """ find local maxima that has a magnitude no smaller than the frac * global maxima. 
        Then merge adjacent peaks, where the maximal height and minimal height between the two peaks are within merge_peak_frac of the maximal height.
    """
    lower_bound = y.max() * frac
    maxima_by_x = []
    filtered_maxima = []
    for i in range(2, y.size - 2):
        if (y[i - 1] == y[i] and y[i - 2] < y[i - 1] and y[i] > y[i + 1]) or (y[i - 2] < y[i - 1] and y[i - 1] < y[i] and y[i] > y[i + 1] and y[i + 1] > y[i + 2]):
            # i is a local maxima
            if y[i] >= lower_bound:
                maxima_by_x.append(i)
            else:
                filtered_maxima.append(i)
    maxima_by_x = np.array(maxima_by_x)
    filtered_maxima = np.array(filtered_maxima)
    n_max = maxima_by_x.size

    curr_peak = 0
    merged_peaks = []
    for i in range(n_max - 1):
        min_value = y[maxima_by_x[i]+1:maxima_by_x[i + 1]].min()
        max_value = max(y[maxima_by_x[i]], y[maxima_by_x[i + 1]])
        if (max_value - min_value) / max_value > merge_peak_frac: # do not merge i + 1
            merged_peaks.append(maxima_by_x[curr_peak])
            curr_peak = i + 1
        else:
            if y[maxima_by_x[i + 1]] > y[maxima_by_x[curr_peak]]:
                curr_peak = i + 1
    merged_peaks.append(maxima_by_x[curr_peak])
    merged_peaks = np.array(merged_peaks)
    maxima = merged_peaks[np.argsort(y[merged_peaks])[::-1]]

    return maxima, maxima_by_x, filtered_maxima

def _find_pos_curv(curv, start, dir, thre = 0.06):
    RANGE = range(start, curv.size) if dir == '+' else range(start, 0, -1)
    assert (RANGE.stop - RANGE.start) * RANGE.step > 0
    for pos in RANGE:
        if curv[pos] > thre:
            break
    return pos

def _find_curv_minima_at_peak(curv, peak_pos):
    start = peak_pos
    while start > 1 and curv[start] < 0.0:
        start -= 1
    start += 1
    end = peak_pos
    while end < curv.size - 2 and curv[end] < 0.0:
        end += 1
    return curv[start:end].min()

def _find_curv_local_minima(curv, peak_curv_value, filtered_maxima, start, dir, rel_thre = 0.4, minima_dir_thre = -0.25):
    """ Find a negative curvature value that is a local minima or a filtered local maxima with respect to density value.
        dir represents the direction of search, choosing from '+' or '-'.
        Beside being a local minima, the value must also satisfy the rel_thre requirement.
        rel_thre requires that the curvature value must smaller than rel_thre fraction of the max of minimal curvature value of the peak and the minimal curvature value since start at direction dir.
    """
    if dir == '+':
        pos_from = max(start, 2)
        pos_to = curv.size - 2
        tmp_arr = filtered_maxima[filtered_maxima > start]
        if tmp_arr.size > 0:
            lmax = tmp_arr.min()
            pos_to = _find_pos_curv(curv, lmax-1, '-') + 1
        assert pos_from < pos_to
        minima_with_dir = curv[pos_from:pos_to].min()
        if minima_with_dir >= minima_dir_thre:
            # No other local minima
            return pos_to # return right end
        thre = max(peak_curv_value, minima_with_dir) * rel_thre
        assert thre < 0.0
        for pos in range(pos_from, pos_to):
            if curv[pos] < thre and curv[pos - 1] > curv[pos] and curv[pos] < curv[pos + 1]:
                return pos
        assert False
    else:
        assert dir == '-'
        pos_from = min(start, curv.size - 3)
        pos_to = 1
        tmp_arr = filtered_maxima[filtered_maxima < start]
        if tmp_arr.size > 0:
            lmax = tmp_arr.max()
            pos_to = _find_pos_curv(curv, lmax+1, '+') - 1
        assert pos_from > pos_to
        minima_with_dir = curv[pos_to+1:pos_from+1].min()
        if minima_with_dir >= minima_dir_thre:
            return pos_to
        thre = max(peak_curv_value, minima_with_dir) * rel_thre
        assert thre < 0.0
        for pos in range(pos_from, pos_to, -1):
            if curv[pos] < thre and curv[pos - 1] > curv[pos] and curv[pos] < curv[pos + 1]:
                return pos
        assert False

def _plot_hist(obs_scores, sim_scores, threshold, sim_x, sim_y, curv, nbin = 100, fig_size = (8,6), dpi = 300):
    """ Plot histogram of doublet scores for observed cells and simulated doublets
        (A) top left: histogram of observed cells;
        (B) top right: histogram of simulated doublets;
        (C) bottom left: KDE of simulated doublets scores
        (D) bottom right: KDE of simulated doublets in log scale
    """
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize = fig_size, dpi = dpi)

    x = np.linspace(0, 1, nbin)
    ax = axes[0, 0]
    ax.hist(obs_scores, x, color="gray", linewidth=0, density=True)
    ax.set_yscale("log")
    ax.axvline(x = threshold, ls = "--", c = "k", linewidth=1)
    ax.set_title('Observed cells')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

    ax = axes[0, 1]
    ax.hist(sim_scores, x, color="gray", linewidth=0, density=True)
    ax.set_yscale("log")
    ax.axvline(x = threshold, ls = "--", c = "k", linewidth=1)
    ax.set_title('Simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

    # ax = axes[1, 0]
    # from scipy.stats import gaussian_kde
    # kde = gaussian_kde(sim_scores)
    # y = kde(x)
    # ax.plot(x, y, '-', c='k', lw = 1)
    # ax.set_ylim(bottom = 0.0)
    # ax.set_title('KDE of simulated doublets')
    # ax.set_xlabel('Doublet score')
    # ax.set_ylabel('Density')

    ax = axes[1, 0]
    ax.plot(sim_x, sim_y, '-', c='k', lw = 1)
    ax.set_ylim(bottom = 0.0)
    ax.axvline(x = np.log(threshold), ls = "--", c="k", lw=1)
    ax.set_title('KDE of simulated doublets')
    ax.set_xlabel('Log doublet score')
    ax.set_ylabel('Density')

    ax = axes[1, 1]
    ax.plot(sim_x, curv, '-', c='k', lw = 1)
    ax.axvline(x = np.log(threshold), ls = "--", c="k", lw=1)
    ax.set_title('Curvature of simulated doublets')
    ax.set_xlabel('Log doublet score')
    ax.set_ylabel('Curvature')

    fig.tight_layout()
    return fig

def _calc_expected_doublet_rate(ncells):
    """ Calculate expected doublet rate based number of cells using 10x Genomics' doublet table [https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-].
        Poisson lambda estimated from table is lambda = 0.00785
    """
    ncell_base = 500.0
    lmd_base = 0.00785

    lmd = lmd_base * (ncells / ncell_base)
    expected_rate = (1.0 - (1.0 + lmd) * np.exp(-lmd)) / (1.0 - np.exp(-lmd))

    return expected_rate


@timer(logger=logger)
def _run_scrublet(
    data: Union[MultimodalData, UnimodalData],
    name: Optional[str] = '',
    expected_doublet_rate: Optional[float] = None,
    sim_doublet_ratio: Optional[float] = 2.0,
    n_prin_comps: Optional[int] = 30,
    robust: Optional[bool] = False,
    k: Optional[int] = None,
    n_jobs: Optional[int] = -1,
    random_state: Optional[int] = 0,
    plot_hist: Optional[bool] = True
) -> Union[None, "Figure"]:
    """Calculate doublet scores using Scrublet-like [Wolock18]_ strategy for the current data.X; determine a right threshold using Gaussian Mixture model.
       This function should be called after highly_variable_gene selection.

    Parameters
    -----------
    data: ``Union[MultimodalData, UnimodalData]`` object.
        Annotated data matrix with rows for cells and columns for genes. Data must be low quality cell and gene filtered and log-transformed. Assume 'raw.X' stores the raw count matrix.

    name: ``str``, optional, default: ``''``
        Name of the sample.

    expected_doublet_rate: ``float``, optional, default: ``None``
        The expected doublet rate for the experiment. By default, calculate the expected rate based on number of cells from the 10x multiplet rate table

    sim_doublet_ratio: ``float``, optional, default: ``2.0``
        The ratio between synthetic doublets and observed cells.

    n_prin_comps: ``int``, optional, default: ``30``
        Number of principal components.

    robust: ``bool``, optional, default: ``False``.
        If true, use 'arpack' instead of 'randomized' for large matrices (i.e. max(X.shape) > 500 and n_components < 0.8 * min(X.shape))

    k: ``int``, optional, default: ``None``
        Number of observed cell neighbors. If None, k = round(0.5 * sqrt(number of observed cells)). Total neighbors k_adj = round(k * (1.0 + sim_doublet_ratio)).

    n_job: ``int``, optional, default: ``-``
        Number of threads to use. If ``-1``, use all available threads.

    random_state: ``int``, optional, default: ``0``
        Random state for doublet simulation, PCA and approximate nearest neighbor search.

    plot_hist: ``bool``, optional, default: ``True``
        If True, plot diagnostic histogram.

    Returns
    --------
    ``None`` or a ``matplotlib Figure object`` if

    Update ``data.obs``:
        * ``data.obs['doublet_score']``: The calculated doublet scores on cells.
        * ``data.obs['pred_dbl']``: Predicted doublets as True.

    Update ``data.uns``:
        * ``data.uns['doublet_threshold']``: Inferred doublet threshold; any score > threshold is identified as a neotypic doublet.

    Examples
    --------
    >>> pg.run_scrublet(data)
    """
    from pegasus.tools import calculate_nearest_neighbors
    from pegasus.cylib.fast_utils import simulate_doublets
    from sklearn.decomposition import PCA
    from scipy.stats import gaussian_kde

    if "highly_variable_features" not in data.var:
        raise ValueError("_run_scrublet must be run after highly_variable_features is called!")

    r = sim_doublet_ratio
    if expected_doublet_rate is None:
        expected_doublet_rate = _calc_expected_doublet_rate(data.shape[0])
    rho = expected_doublet_rate

    # subset the raw count matrix
    rawX = data.get_matrix("raw.X")
    obs_umis = rawX.sum(axis = 1, dtype = np.int32).A1
    rawX = rawX[:, data.var["highly_variable_features"].values]
    # Simulate synthetic doublets
    sim_rawX, pair_idx = simulate_doublets(rawX, r, random_state)
    sim_umis = obs_umis[pair_idx].sum(axis = 1, dtype = np.int32)

    # standardize and calculate PCA for rawX
    obsX = rawX.astype(np.float32).toarray()
    obsX /= obs_umis.reshape(-1, 1) # normalize each cell

    m1 = obsX.mean(axis = 0) # calculate mean and std
    psum = np.multiply(obsX, obsX).sum(axis=0)
    std = ((psum - obsX.shape[0] * (m1 ** 2)) / (obsX.shape[0] - 1.0)) ** 0.5
    std[std == 0] = 1

    obsX -= m1 # standardize
    obsX /= std

    svd_solver = "auto" if not robust else ("arpack" if max(obsX.shape) > 500 and n_prin_comps < 0.8 * min(obsX.shape) else "full") # PCA
    pca = PCA(n_components=n_prin_comps, random_state=random_state, svd_solver=svd_solver)
    obs_pca = pca.fit_transform(obsX)

    # standardize and calculate PCA for sim_rawX
    simX = sim_rawX.astype(np.float32).toarray()
    simX /= sim_umis.reshape(-1, 1) # normalize each cell

    simX -= m1 # standardize
    simX /= std

    sim_pca = pca.transform(simX) # transform to PC coordinates

    # concatenate observed and simulated data
    pc_coords = np.vstack((obs_pca, sim_pca))
    is_doublet = np.repeat(np.array([0, 1], dtype = np.int32), [obsX.shape[0], simX.shape[0]])

    # Calculate k nearest neighbors
    if k is None:
        k = int(round(0.5 * np.sqrt(obsX.shape[0])))
    k_adj = int(round(k * (1.0 + r)))
    indices, _ = calculate_nearest_neighbors(pc_coords, K = k_adj + 1, n_jobs = n_jobs)

    # Calculate scrublet-like doublet score
    k_d = is_doublet[indices].sum(axis = 1)
    q = (k_d + 1.0) / (k_adj + 2.0) # Equation 5
    doublet_scores = (q * rho / r) / ((1.0 - rho) - q * (1.0 - rho - rho / r)) # Equation 4
    obs_scores = doublet_scores[0:obsX.shape[0]]
    sim_scores = doublet_scores[obsX.shape[0]:]

    # Determine a scrublet score threshold
    # log transformed
    sim_scores_log = np.log(sim_scores)

    # Estimate KDE
    min_score = sim_scores_log.min()
    max_score = sim_scores_log.max()
    min_gap = np.diff(np.unique(np.sort(sim_scores_log))).min()
    from math import ceil
    n_gap = max(int(ceil((max_score - min_score) / min_gap)), 200) # minimum is 200
    gap = (max_score - min_score) / n_gap

    n_ext = 5
    min_score -= gap * n_ext
    max_score += gap * n_ext
    x = np.linspace(min_score, max_score, n_gap + 1 + n_ext * 2) # generate x coordinates
    kde = gaussian_kde(sim_scores_log)
    y = kde(x)

    # Find local maxima
    maxima, maxima_by_x, filtered_maxima = _find_local_maxima(y)
    assert maxima.size > 0
    curv = _calc_vec_f(_curvature, x.size, y, gap) # calculate curvature

    if maxima.size >= 2:
        if maxima[0] < maxima[1]:
            start = maxima[0]
            end = maxima[1]
        else:
            start = maxima[1]
            end = maxima[0]
        pos = y[start+1:end].argmin() + (start+1)
    else:
        frac_right_thre = 0.42
        frac_left_thre = 0.4

        pos = -1
        for i in range(maxima_by_x.size):
            frac_right = (sim_scores_log > x[maxima_by_x[i]]).sum() / sim_scores.size
            if frac_right < frac_right_thre: # peak might represent a doublet peak, try to find a cutoff at the left side
                if i == 0:
                    peak_curv_value = _find_curv_minima_at_peak(curv, maxima_by_x[i])
                    end = _find_pos_curv(curv, maxima_by_x[i]-1, '-')
                    start = _find_pos_curv(curv, _find_curv_local_minima(curv, peak_curv_value, filtered_maxima, end-1, '-')+1, '+')
                    assert start <= end
                    pos = curv[start:end+1].argmax() + start
                else:
                    pos = y[maxima_by_x[i-1]+1:maxima_by_x[i]].argmin() + (maxima_by_x[i-1]+1)

                frac_left = (sim_scores_log < x[pos]).sum() / sim_scores.size    
                if frac_left < frac_left_thre:
                    pos = maxima_by_x[i]

                break

        if pos < 0:
            # peak represents singlet, find a cutoff at the right side
            peak_curv_value = _find_curv_minima_at_peak(curv, maxima_by_x[-1])
            start = _find_pos_curv(curv, maxima_by_x[-1]+1, '+')
            end = _find_pos_curv(curv, _find_curv_local_minima(curv, peak_curv_value, filtered_maxima, start+1, '+')-1, '-')
            assert start <= end
            pos = curv[start:end+1].argmax() + start

    threshold = np.exp(x[pos])

    data.obs["doublet_score"] = obs_scores.astype(np.float32)
    data.obs["pred_dbl"] = obs_scores > threshold
    data.uns["doublet_threshold"] = float(threshold)

    logger.info(f"Sample {name}: doublet threshold = {threshold:.4f}; total cells = {data.shape[0]}; neotypic doublet rate = {data.obs['pred_dbl'].sum() / data.shape[0]:.2%}")
    fig = None
    if plot_hist:
        fig = _plot_hist(obs_scores, sim_scores, threshold, x, y, curv)
    return fig


def _identify_doublets_fisher(cluster_labels: Union[pd.Categorical, List[int]], pred_dbl: List[bool], alpha: float = 0.05) -> pd.DataFrame:
    df = pd.crosstab(cluster_labels, pred_dbl)

    ndbl = df[True].sum()
    a = df[True].values.astype(np.int32)
    b = df[False].values.astype(np.int32)
    c = ndbl - a
    d = (pred_dbl.size - ndbl) - b

    avg_dblr = ndbl / pred_dbl.size
    freqs = a / (a + b)

    from pegasus.cylib.cfisher import fisher_exact
    from statsmodels.stats.multitest import fdrcorrection as fdr
    _, pvals = fisher_exact(a, b, c, d)
    passed, qvals = fdr(pvals, alpha = alpha)

    posvec = np.where(passed)[0][freqs[passed] > avg_dblr]

    result = pd.DataFrame({'cluster': df.index[posvec], 'percentage': freqs[posvec] * 100.0, 'pval': pvals[posvec], 'qval': qvals[posvec]})
    result.sort_values('percentage', ascending = False, inplace = True)
    result.reset_index(drop=True, inplace=True)

    return result



@timer(logger=logger)
def infer_doublets(
    data: MultimodalData,
    channel_attr: Optional[str] = None,
    clust_attr: Optional[str] = None,
    min_cell: Optional[int] = 100,
    expected_doublet_rate: Optional[float] = None,
    sim_doublet_ratio: Optional[float] = 2.0,
    n_prin_comps: Optional[int] = 30,
    robust: Optional[bool] = False,
    k: Optional[int] = None,
    n_jobs: Optional[int] = -1,
    alpha: Optional[float] = 0.05,
    random_state: Optional[int] = 0,
    plot_hist: Optional[str] = "dbl",
) -> None:
    """Infer doublets using a Scrublet-like strategy. [Li20-2]_

    This function must be called after clustering. 

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    channel_attr: ``str``, optional, default: None
        Attribute indicating sample channels. If set, calculate scrublet-like doublet scores per channel.

    clust_attr: ``str``, optional, default: None
        Attribute indicating cluster labels. If set, estimate proportion of doublets in each cluster and statistical significance.

    min_cell: ``int``, optional, default: 100
        Minimum number of cells per sample to calculate doublet scores. For samples having less than 'min_cell' cells, doublet score calculation will be skipped.

    expected_doublet_rate: ``float``, optional, default: ``None``
        The expected doublet rate for the experiment. By default, calculate the expected rate based on number of cells from the 10x multiplet rate table

    sim_doublet_ratio: ``float``, optional, default: ``2.0``
        The ratio between synthetic doublets and observed cells.

    n_prin_comps: ``int``, optional, default: ``30``
        Number of principal components.

    robust: ``bool``, optional, default: ``False``.
        If true, use 'arpack' instead of 'randomized' for large matrices (i.e. max(X.shape) > 500 and n_components < 0.8 * min(X.shape))

    k: ``int``, optional, default: ``None``
        Number of observed cell neighbors. If None, k = round(0.5 * sqrt(number of observed cells)). Total neighbors k_adj = round(k * (1.0 + sim_doublet_ratio)).

    n_job: ``int``, optional, default: ``-``
        Number of threads to use. If ``-1``, use all available threads.

    alpha: ``float``, optional, default: ``0.05``
        FDR significant level for cluster-level fisher exact test.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    plot_hist: ``str``, optional, default: ``dbl``
        If not None, plot diagnostic histograms using ``plot_hist`` as the prefix. If `channel_attr` is None, ``plot_hist.png`` is generated; Otherwise, ``plot_hist.channel_name.png`` files are generated.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs['pred_dbl_type']``: Predicted singlet/doublet types.

        * ``data.uns['pred_dbl_cluster']``: Only generated if 'clust_attr' is not None. This is a dataframe with two columns, 'Cluster' and 'Qval'. Only clusters with significantly more doublets than expected will be recorded here.

    Examples
    --------
    >>> pg.infer_doublets(data, channel_attr = 'Channel', clust_attr = 'Annotation')
    """
    assert data.get_modality() == "rna"
    try:
        rawX = data.get_matrix("raw.X")
    except ValueError:
        raise ValueError("Cannot detect the raw count matrix raw.X; stop inferring doublets!")

    if_plot = plot_hist is not None

    if channel_attr is None:
        if data.shape[0] >= min_cell:
            fig = _run_scrublet(data, expected_doublet_rate = expected_doublet_rate, sim_doublet_ratio = sim_doublet_ratio, \
                                n_prin_comps = n_prin_comps, robust = robust, k = k, n_jobs = n_jobs, random_state = random_state, \
                                plot_hist = if_plot)
            if if_plot:
                fig.savefig(f"{plot_hist}.png")
        else:
            logger.warning(f"Data has {data.shape[0]} < {min_cell} cells and thus doublet score calculation is skipped!")
            data.obs["doublet_score"] = 0.0
            data.obs["pred_dbl"] = False
    else:
        from pandas.api.types import is_categorical_dtype
        from pegasus.tools import identify_robust_genes, log_norm, highly_variable_features

        assert is_categorical_dtype(data.obs[channel_attr])
        genome = data.get_genome()
        modality = data.get_modality()
        channels = data.obs[channel_attr].cat.categories

        dbl_score = np.zeros(data.shape[0], dtype = np.float32)
        pred_dbl = np.zeros(data.shape[0], dtype = np.bool_)
        thresholds = {}
        for channel in channels:
            # Generate a new unidata object for the channel
            idx = np.where(data.obs[channel_attr] == channel)[0]
            if idx.size >= min_cell:
                unidata = UnimodalData({"barcodekey": data.obs_names[idx]}, 
                                       {"featurekey": data.var_names},
                                       {"X": rawX[idx]},
                                       {"genome": genome, "modality": modality})
                # Identify robust genes, count and log normalized and select top 2,000 highly variable features
                identify_robust_genes(unidata)
                log_norm(unidata)
                highly_variable_features(unidata)
                # Run _run_scrublet
                fig = _run_scrublet(unidata, name = channel, expected_doublet_rate = expected_doublet_rate, sim_doublet_ratio = sim_doublet_ratio, \
                                    n_prin_comps = n_prin_comps, robust = robust, k = k, n_jobs = n_jobs, random_state = random_state, \
                                    plot_hist = if_plot)
                if if_plot:
                    fig.savefig(f"{plot_hist}.{channel}.png")

                dbl_score[idx] = unidata.obs["doublet_score"].values
                pred_dbl[idx] = unidata.obs["pred_dbl"].values
                thresholds[channel] = unidata.uns["doublet_threshold"]
            else:
                logger.warning(f"Channel {channel} has {idx.size} < {min_cell} cells and thus doublet score calculation is skipped!")

        data.obs["doublet_score"] = dbl_score
        data.obs["pred_dbl"] = pred_dbl
        data.uns["doublet_thresholds"] = thresholds

    if clust_attr is not None:
        data.uns["pred_dbl_cluster"] = _identify_doublets_fisher(data.obs[clust_attr].values, data.obs["pred_dbl"].values, alpha = alpha)

    logger.info('Doublets are predicted!')



def mark_doublets(
    data: MultimodalData,
    demux_attr: Optional[str] = 'demux_type',
    dbl_clusts: Optional[str] = None,
) -> None:
    """Convert doublet prediction into doublet annotations that Pegasus can recognize. In addition, clusters in dbl_clusts will be marked as doublets.

    Must run ``infer_doublets`` first.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    demux_attr: ``str``, optional, default: ``demux_type``
        Attribute indicating singlets/doublets that Pegasus can recognize. Currently this is 'demux_type', which is also used for hashing.

    dbl_clusts: ``str``, optional, default: None
        Indicate which clusters should be marked as all doublets. It takes the format of 'clust:value1,value2,...', where 'clust' refers to the cluster attribute.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[demux_attr]``: Singlet/doublet annotation.

    Examples
    --------
    >>> pg.mark_doublets(data, dbl_clusts='Annotation:B/T doublets')
    """
    codes = data.obs["pred_dbl"].values.astype(np.int32)
    if dbl_clusts is not None:
        cluster, value_str = dbl_clusts.split(':')
        idx = np.isin(data.obs[cluster], value_str.split(','))
        codes[idx] = 1
    data.obs[demux_attr] = pd.Categorical.from_codes(codes, categories = ["singlet", "doublet"])
