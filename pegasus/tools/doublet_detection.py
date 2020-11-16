import time
import numpy as np
import pandas as pd

from typing import List, Optional, Union

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

def _find_local_maxima(y, fc = 6.6): # find local maxima that has a magnitude no smaller than the global maxima / fc.
    lower_bound = 0.0
    posvec = np.argsort(y[2:y.size - 2])[::-1] + 2

    maxima = []
    for i in posvec:
        if y[i] < lower_bound:
            break
        if (y[i - 1] == y[i] and y[i - 2] < y[i - 1] and y[i] > y[i + 1]) or (y[i - 2] < y[i - 1] and y[i - 1] < y[i] and y[i] > y[i + 1] and y[i + 1] > y[i + 2]):
            if lower_bound == 0.0: # y[i] is the maximum value
                lower_bound = y[i] / fc
            maxima.append(i)

    return np.array(maxima)

def _find_pos_curv(curv, start, dir):
    RANGE = range(start, curv.size) if dir == '+' else range(start, 0, -1)
    for pos in RANGE:
        if curv[pos] > 0.0:
            break
    return pos

def _find_curv_local_minima(curv, start, dir, thre = -1.0):
    """ find a negative curvature value that is a local minima and no larger than thre
        from "max(start, 2)" to "curv.size - 2", making sure curvature is not nan
        dir '+' or '-'
    """
    RANGE = range(max(start, 2), curv.size - 2) if dir == '+' else range(min(start, curv.size - 3), 1, -1)
    for pos in RANGE:
        if curv[pos] <= thre and curv[pos - 1] > curv[pos] and curv[pos] < curv[pos + 1]:
            break
    return pos

def _plot_hist(obs_scores, sim_scores, threshold, sim_x, sim_y, nbin = 100, fig_size = (8,6), dpi = 300):
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

    ax = axes[1, 0]
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(sim_scores)
    y = kde(x)
    ax.plot(x, y, '-', c='k', lw = 1)
    ax.set_ylim(bottom = 0.0)
    ax.set_title('KDE of simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

    ax = axes[1, 1]
    ax.plot(sim_x, sim_y, '-', c='k', lw = 1)
    ax.set_ylim(bottom = 0.0)
    ax.axvline(x = np.log(threshold), ls = "--", c="k", lw=1)
    ax.set_title('KDE of simulated doublets (log scale)')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

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
    from sklearn.mixture import GaussianMixture
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
    sslog_reshaped = sim_scores_log.reshape(-1, 1)

    # Gaussian Mixture
    gm = GaussianMixture(n_components = 3, random_state = random_state)
    gm.fit(sslog_reshaped)
    preds = gm.predict(sslog_reshaped)
    dbl_code = np.argmax(gm.means_.ravel())
    gm_thre = (sim_scores_log[preds != dbl_code].max() + sim_scores_log[preds == dbl_code].min()) / 2.0 # boundary between singlets and doublets as inferred by the GM model

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
    maxima = _find_local_maxima(y)
    assert maxima.size > 0

    slt_maxima = maxima[x[maxima] < gm_thre]
    dbl_maxima = maxima[x[maxima] > gm_thre]

    if slt_maxima.size > 0 and dbl_maxima.size > 0:
        # We have peaks in both singlet and doublet portions
        pos = y[slt_maxima.max()+1:dbl_maxima.min()].argmin() + (slt_maxima.max()+1)
    else:
        curv = _calc_vec_f(_curvature, x.size, y, gap)
        if slt_maxima.size > 0:
            # We have only one peak and peak is in the singlet portion
            leftmost_dbl = np.where(x > gm_thre)[0][0]
            start = _find_pos_curv(curv, max(slt_maxima.max(), _find_curv_local_minima(curv, leftmost_dbl-1, '-'))+1, '+')
            end = _find_pos_curv(curv, _find_curv_local_minima(curv, leftmost_dbl, '+')-1, '-')
        else:
            # We have only one peak and peak is in the doublet portion
            start = _find_pos_curv(curv, dbl_maxima.min()+1, '+')
            end = _find_pos_curv(curv, _find_curv_local_minima(curv, start+1, '+')-1, '-')
        assert start <= end
        pos = np.abs(curv[start:end+1]).argmax() + start

    threshold = np.exp(x[pos])

    data.obs["doublet_score"] = obs_scores.astype(np.float32)
    data.obs["pred_dbl"] = obs_scores > threshold
    data.uns["doublet_threshold"] = float(threshold)

    fig = None
    if plot_hist:
        fig = _plot_hist(obs_scores, sim_scores, threshold, x, y)
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
        fig = _run_scrublet(data, expected_doublet_rate = expected_doublet_rate, sim_doublet_ratio = sim_doublet_ratio, \
                            n_prin_comps = n_prin_comps, robust = robust, k = k, n_jobs = n_jobs, random_state = random_state, \
                            plot_hist = if_plot)
        if if_plot:
            fig.savefig(f"{plot_hist}.png")
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
            unidata = UnimodalData({"barcodekey": data.obs_names[idx]},
                                   {"featurekey": data.var_names},
                                   {"X": rawX[idx]},
                                   {"genome": genome, "modality": modality})
            # Identify robust genes, count and log normalized and select top 2,000 highly variable features
            identify_robust_genes(unidata)
            log_norm(unidata)
            highly_variable_features(unidata)
            # Run _run_scrublet
            fig = _run_scrublet(unidata, expected_doublet_rate = expected_doublet_rate, sim_doublet_ratio = sim_doublet_ratio, \
                                n_prin_comps = n_prin_comps, robust = robust, k = k, n_jobs = n_jobs, random_state = random_state, \
                                plot_hist = if_plot)
            if if_plot:
                fig.savefig(f"{plot_hist}.{channel}.png")

            dbl_score[idx] = unidata.obs["doublet_score"].values
            pred_dbl[idx] = unidata.obs["pred_dbl"].values
            thresholds[channel] = unidata.uns["doublet_threshold"]

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
