import time
import numpy as np
import pandas as pd

from scipy.sparse import issparse, csr_matrix
from scipy.signal import find_peaks
from typing import List, Optional, Union, Tuple, Dict
from matplotlib.figure import Figure

import logging
logger = logging.getLogger(__name__)

from threadpoolctl import threadpool_limits

from pegasusio import UnimodalData, MultimodalData
from pegasusio import timer
from pegasus.tools import eff_n_jobs



def _calc_expected_doublet_rate(ncells):
    """ Calculate expected doublet rate based number of cells using 10x Genomics' doublet table [https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-].
        Poisson lambda estimated from table is lambda = 0.00785
    """
    ncell_base = 500.0
    lmd_base = 0.00785

    lmd = lmd_base * (ncells / ncell_base)
    expected_rate = (1.0 - (1.0 + lmd) * np.exp(-lmd)) / (1.0 - np.exp(-lmd))

    return expected_rate


def _simulate_doublets(X: Union[csr_matrix, np.ndarray], sim_doublet_ratio: float, random_state: int = 0) -> Tuple[Union[csr_matrix, np.ndarray], np.ndarray]:
    """ Function to simulate doublets by randomly picking two observed cells and merge them togethers
    """
    # simulate doublet indices
    np.random.seed(random_state)
    n_sim = int(X.shape[0] * sim_doublet_ratio)
    doublet_indices = np.random.randint(0, X.shape[0], size=(n_sim, 2), dtype = np.int32)

    results = None
    if issparse(X):
        data = X.data
        if data.dtype != np.int32:
            data = data.astype(np.int32)
        from pegasus.cylib.fast_utils import simulate_doublets_sparse
        results = csr_matrix(simulate_doublets_sparse(n_sim, X.shape[1], data, X.indices, X.indptr, doublet_indices), shape = (n_sim, X.shape[1]), copy = False)
    else:
        data = X
        if data.dtype != np.int32:
            data = data.astype(np.int32)
        from pegasus.cylib.fast_utils import simulate_doublets_dense
        results = simulate_doublets_dense(n_sim, X.shape[1], data, doublet_indices)

    return results, doublet_indices


def _gen_sim_data(
    data: MultimodalData,
    raw_mat_key: str,
    sim_doublet_ratio: float,
    random_state: int,
    anno_attr: Optional[str] = None,
) -> MultimodalData:

    rawX = data.get_matrix(raw_mat_key)
    simX, pair_idx = _simulate_doublets(rawX, sim_doublet_ratio, random_state)

    barcodes = data.obs_names.values.astype(str)
    barcodekey = np.char.add(np.char.add(barcodes[pair_idx[:, 0]], '|'), barcodes[pair_idx[:, 1]])

    import pegasusio as pio

    unidat = pio.UnimodalData({'barcodekey': barcodekey},
             {'featurekey': data.var_names},
             {'counts': simX},
             {'modality': data.uns['modality'], 'genome': data.uns['genome']})

    if 'featureid' in data.var:
        unidat.var['featureid'] = data.var['featureid']
    assert ('robust' in data.var) and ('highly_variable_features' in data.var)
    unidat.var['robust'] = data.var['robust']
    unidat.var['highly_variable_features'] = data.var['highly_variable_features']

    if anno_attr is not None:
        assert anno_attr in data.obs
        annos = data.obs[anno_attr].values.astype(str)
        arr = np.sort(np.column_stack((annos[pair_idx[:, 0]], annos[pair_idx[:, 1]])))
        unidat.obs['anno'] = arr[:, 0]
        idx = arr[:, 0] != arr[:, 1]
        unidat.obs.loc[idx, 'anno'] = np.char.add(np.char.add(arr[idx, 0], '|'), arr[idx, 1])
        unidat.obs['anno'] = pd.Categorical(unidat.obs['anno'])

    sim_data = pio.MultimodalData(unidat)

    from pegasus.tools import qc_metrics, log_norm

    qc_metrics(sim_data)
    assert 'norm_count' in data.uns
    log_norm(sim_data, norm_count=data.uns['norm_count'])

    return sim_data

def _concat_and_process(datslt, datsim, K, r, rho, n_jobs, calc_umap=False, anno_attr=None):
    from scipy.sparse import vstack
    import pegasusio as pio
    from pegasus.tools import pc_transform, neighbors, umap

    assert ('robust' in datslt.var) and ('highly_variable_features' in datslt.var)
    assert ('n_genes' in datslt.obs) and ('n_counts' in datslt.obs) and ('n_genes' in datsim.obs) and ('n_counts' in datsim.obs)

    unidat = pio.UnimodalData({'barcodekey': np.concatenate((datslt.obs_names, datsim.obs_names)),
                               'type': np.repeat(('real', 'sim'), (datslt.shape[0], datsim.shape[0])),
                               'n_genes': np.concatenate((datslt.obs['n_genes'].values, datsim.obs['n_genes'].values)),
                               'n_counts': np.concatenate((datslt.obs['n_counts'].values, datsim.obs['n_counts'].values)),
                              },
             {'featurekey': datslt.var_names,
              'robust': datslt.var['robust'],
              'highly_variable_features': datslt.var['highly_variable_features'],
             },
             {'counts': vstack([datslt.get_matrix('counts'), datsim.get_matrix('counts')]),
              'counts.log_norm': vstack([datslt.get_matrix('counts.log_norm'), datsim.get_matrix('counts.log_norm')]),
             },
             {'modality': datslt.uns['modality'], 'genome': datslt.uns['genome'], 'pca_ncomps': datslt.uns['pca_ncomps'] if "pca_ncomps" in datslt.uns else datslt.obsm["X_pca"].shape[1]},
             cur_matrix='counts.log_norm')

    if 'featureid' in datslt.var:
        unidat.var['featureid'] = datslt.var['featureid']

    if (anno_attr is not None) and (anno_attr in datslt.obs):
        unidat.obs['anno'] = np.concatenate((datslt.obs[anno_attr], np.repeat('sim', datsim.shape[0])))
        if 'anno' in datsim.obs:
            unidat.obs['anno_all'] = np.concatenate((datslt.obs[anno_attr], np.char.add('!', datsim.obs['anno'].values.astype(str))))
            unidat.obs['subtype'] = unidat.obs['type'].copy()
            unidat.obs.loc[unidat.obs['anno_all'].map(lambda x: x.find('|') > 0), 'subtype'] = 'sim_neo'

    unidat.obsm['X_pca'] = np.vstack((datslt.obsm['X_pca'], pc_transform(datslt, datsim.get_matrix('counts'))))
    dat_concat = pio.MultimodalData(unidat)

    key = 'pca'
    n_comps = dat_concat.uns['pca_ncomps']
    neighbors(dat_concat, K=K, rep=key, n_comps=n_comps, n_jobs=n_jobs)
    if calc_umap:
        umap(dat_concat, rep=key, rep_ncomps=n_comps, n_jobs=n_jobs)

    knn_key = f'{key}_knn_indices'
    k = dat_concat.obsm[knn_key].shape[1]

    dat_concat.obs['q'] = ((dat_concat.obs['type'].values[dat_concat.obsm[knn_key]] == 'sim').sum(axis=1) + 1.0) / (k + 2.0) # Equation 5
    dat_concat.obs['doublet_score'] = (dat_concat.obs['q'] * rho / r) / ((1.0 - rho) - dat_concat.obs['q'] * (1.0 - rho - rho / r)) # Equation 4

    return dat_concat


def _calc_expected_emb_rate_in_sim(obs_pca, rho, n_jobs, random_state):
    """ What's the expected embedding doublet rate if assume Kmeans with 5 cluster as the ground truth cell type? This will give us some guidance for determining doublet cutoff
    """
    from sklearn.cluster import KMeans

    with threadpool_limits(limits = n_jobs):
        kmeans = KMeans(n_clusters = 5, random_state = random_state, n_init='auto').fit(obs_pca)

    # calculate in simulated distribution, expected percentage of embedded doublets
    _, freqs = np.unique(kmeans.labels_, return_counts = True)
    freqs = np.array(freqs) / sum(freqs)
    d_emb = (((1.0 - rho) * freqs + rho * (freqs ** 2)) ** 2).sum()
    d_neo = 1.0 - d_emb

    return kmeans.labels_, d_emb, d_neo



def _kde_smooth(scores, bw_method='silverman'):
    from scipy.stats import gaussian_kde
    from math import ceil

    # Estimate KDE
    min_score = scores.min()
    max_score = scores.max()
    min_gap = np.diff(np.unique(np.sort(scores))).min()
    n_gap = max(int(ceil((max_score - min_score) / min_gap)), 200) # minimum is 200
    gap = (max_score - min_score) / n_gap

    n_ext = 5
    min_score -= gap * n_ext
    max_score += gap * n_ext
    x = np.linspace(min_score, max_score, n_gap + 1 + n_ext * 2) # generate x coordinates
    kde = gaussian_kde(scores, bw_method=bw_method)
    y = kde(x)

    return x, y, gap


def _find_local_maxima(y: List[float], frac: float = 0.25, merge_peak_frac: float = 0.06) -> Tuple[List[int], List[int], List[int]]:
    """ find local maxima that has a magnitude larger than the frac * global maxima.
        Then merge adjacent peaks, where the maximal height and minimal height between the two peaks are within merge_peak_frac of the maximal height.

        RETURNS:
          maxima: merged peak coordinates sorted y in descending order
          peak_groups: map of merged peak coordinates to their corresponding unmerged peak coordinates, which all pass the filter
          filtered_maxima: filtered peak coordinates
    """
    lower_bound = y.max() * frac
    all_peaks, _ = find_peaks(y)
    is_above = y[all_peaks] > lower_bound
    maxima_by_x = all_peaks[is_above]
    filtered_maxima = all_peaks[~is_above]
    n_max = maxima_by_x.size

    peak_groups = dict()

    curr_peak = 0
    merged_peaks = []
    tmp_peaks = [maxima_by_x[0]]
    for i in range(n_max - 1):
        min_value = y[maxima_by_x[i]+1:maxima_by_x[i + 1]].min()
        max_value = max(y[maxima_by_x[i]], y[maxima_by_x[i + 1]])
        if (max_value - min_value) / max_value > merge_peak_frac: # do not merge i + 1
            merged_peaks.append(maxima_by_x[curr_peak])
            peak_groups[maxima_by_x[curr_peak]] = np.array(tmp_peaks, dtype=int)
            tmp_peaks = [maxima_by_x[i + 1]]
            curr_peak = i + 1
        else:
            tmp_peaks.append(maxima_by_x[i + 1])
            if y[maxima_by_x[i + 1]] > y[maxima_by_x[curr_peak]]:
                curr_peak = i + 1
    merged_peaks.append(maxima_by_x[curr_peak])
    peak_groups[maxima_by_x[curr_peak]] = np.array(tmp_peaks, dtype=int)
    merged_peaks = np.array(merged_peaks, dtype=int)
    maxima = merged_peaks[np.argsort(y[merged_peaks])[::-1]]

    return maxima, peak_groups, filtered_maxima

def _locate_cutoff_among_peaks_with_guide(x: List[float], y: List[float], maxima: List[float], scores: List[float], d_neo: float) -> int:
    """
        First pick is maxima[0], goal is to find the second largest peak.
        Pick the second largest peak that the cutoff results in a neotypic doublet rate closest to d_neo
    """
    best_delta = 1e100
    best_pos = -1
    for i in range(1, maxima.size):
        if maxima[0] < maxima[i]:
            start = maxima[0]
            end = maxima[i]
        else:
            start = maxima[i]
            end = maxima[0]
        pos = y[start+1:end].argmin() + (start+1)
        d_prac_neo = (scores > x[pos]).sum() / scores.size
        delta = abs(d_prac_neo - d_neo)
        if best_delta > delta:
            best_delta = delta
            best_pos = pos
    return best_pos


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

def _find_pos_curv(curv, start, dir, err_bound = 0.05):
    """ Find a position that has a positive curvature value """
    RANGE = range(start, curv.size) if dir == '+' else range(start, 0, -1)
    assert (RANGE.stop - RANGE.start) * RANGE.step > 0
    for pos in RANGE:
        if curv[pos] > err_bound:
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

def _find_curv_local_minima(curv, peak_curv_value, filtered_maxima, start, rel_thre = 0.1, minima_dir_thre = -0.25):
    """ Find a negative curvature value that is a local minima or a filtered local maxima with respect to density value at the right hand side of start.
        Beside being a local minima, the value must also satisfy the rel_thre requirement.
        rel_thre requires that the curvature value must smaller than rel_thre fraction of the max of minimal curvature value of the peak and the minimal curvature value since start at direction dir.
    """
    pos_from = max(start, 2)
    pos_to = curv.size - 2
    tmp_arr = filtered_maxima[filtered_maxima > start]
    if tmp_arr.size > 0:
        lmax = tmp_arr.min()
        pos_to = _find_pos_curv(curv, lmax-1, '-') + 1
    if pos_from >= pos_to:
        # No local minima within the range
        return pos_to
    minima_with_dir = curv[pos_from:pos_to].min()
    if minima_with_dir >= minima_dir_thre:
        # No local minima lower than minima_dir_thre
        return pos_to # return right end
    thre = min(max(peak_curv_value, minima_with_dir) * rel_thre, minima_dir_thre)
    assert thre < 0.0
    for pos in range(pos_from, pos_to):
        if curv[pos] < thre and curv[pos - 1] > curv[pos] and curv[pos] <= curv[pos + 1]:
            return pos
    return pos_to - 1

def _find_cutoff_right_side(peak_pos: int, peak_groups: Dict[int, List[int]], curv: List[float], filtered_maxima: List[int]) -> int:
    # Peak represents embedded doublets, find a cutoff at the right side
    peak_curv_value = _find_curv_minima_at_peak(curv, peak_pos)
    peak_pos_rlimit = peak_groups[peak_pos].max()
    start = _find_pos_curv(curv, peak_pos_rlimit+1, '+')
    end = _find_pos_curv(curv, _find_curv_local_minima(curv, peak_curv_value, filtered_maxima, start+1)-1, '-')

    return start + (curv[start:end+1].argmax() if start <= end else 0)

def _find_cutoff_left_side(peak_pos: int, x: List[float], curv: List[float], x_theory: float) -> int:
    # Peak represents a doublet peak and thus we need to find a cutoff at the left side
    end = _find_pos_curv(curv, peak_pos-1, '-')
    start = end
    while start > 2 and x[start] >= x_theory:
        start -= 1
    while start > 2 and not (curv[start - 1] > curv[start] and curv[start] < curv[start + 1]):
        start -= 1
    return start + curv[start:end+1].argmax()

def _find_score_threshold(sim_scores, d_neo, threshold_guide, threshold_expected, manual_correction):
    # Gaussian smoothing with bw_method = silverman (a bit more robust)
    x, y, gap = _kde_smooth(sim_scores)

    # Find local maxima
    maxima, peak_groups, filtered_maxima = _find_local_maxima(y)
    assert maxima.size > 0

    # Calculate curvature
    curv = _calc_vec_f(_curvature, x.size, y, gap)

    # Compute cutoff position
    pos = -1
    if maxima.size >= 2:
        pos = _locate_cutoff_among_peaks_with_guide(x, y, maxima, sim_scores, d_neo)
    else:
        frac_right = (sim_scores > x[maxima[0]]).sum() / sim_scores.size
        if x[maxima[0]] < threshold_expected or frac_right > 0.4:
            pos = _find_cutoff_right_side(maxima[0], peak_groups, curv, filtered_maxima)
        else:
            pos = _find_cutoff_left_side(maxima[0], x, curv, threshold_guide)

    threshold = x[pos]

    threshold_auto = None
    if manual_correction is not None:
        threshold_auto = threshold
        if manual_correction == "peak":
            threshold = maxima[0]
        elif manual_correction == "expected":
            threshold = threshold_expected
        elif manual_correction == "right":
            pos = _find_cutoff_right_side(maxima[0], peak_groups, curv, filtered_maxima)
            threshold = x[pos]
        elif manual_correction == "left":
            pos = _find_cutoff_left_side(maxima[0], x, curv, threshold_guide)
            threshold = x[pos]
        else:
            threshold = float(manual_correction)

    return threshold, threshold_auto, x, y, curv


def _calc_bc_sarle(scores):
    """ Sarle's Biomodality Coefficient for finite samples: https://rdrr.io/github/microbiome/microbiome/man/bimodality_sarle.html
    """
    from scipy.stats import skew, kurtosis
    g = skew(scores)
    k = kurtosis(scores)
    n = scores.size
    return (g ** 2 + 1.0) / (k + 3 * (n - 1) ** 2 / (n - 2) / (n - 3))


def _plot_hist(obs_scores, sim_scores, threshold, threshold_guide, threshold_expected, sim_x, sim_y, curv, nbin = 100, fig_size = (8,6), dpi = 300, threshold_auto = None):
    """ Plot histogram of doublet scores for observed cells and simulated doublets
        (A) top left: histogram of observed cells;
        (B) top right: histogram of simulated doublets;
        (C) bottom left: KDE of simulated doublet scores
        (D) bottom right: Curvature of the KDE of simulated doublet scores
    """
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize = fig_size, dpi = dpi)

    x = np.linspace(0, 1, nbin)
    ax = axes[0, 0]
    ax.hist(obs_scores, x, color="gray", linewidth=0, density=True)
    ax.set_yscale("log")
    ax.axvline(x = threshold_expected, ls = "--", c = "b", linewidth=1)
    ax.axvline(x = threshold_guide, ls = "--", c = "r", linewidth=1)
    if threshold_auto is not None:
        ax.axvline(x = threshold_auto, ls = "--", c = "g", linewidth=1)
    ax.axvline(x = threshold, ls = "--", c = "k", linewidth=1)
    ax.set_title('Observed cells')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

    ax = axes[0, 1]
    ax.hist(sim_scores, x, color="gray", linewidth=0, density=True)
    ax.set_yscale("log")
    ax.axvline(x = threshold_expected, ls = "--", c = "b", linewidth=1)
    ax.axvline(x = threshold_guide, ls = "--", c = "r", linewidth=1)
    if threshold_auto is not None:
        ax.axvline(x = threshold_auto, ls = "--", c = "g", linewidth=1)
    ax.axvline(x = threshold, ls = "--", c = "k", linewidth=1)
    ax.set_title('Simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

    ax = axes[1, 0]
    ax.plot(sim_x, sim_y, '-', c='k', lw = 1)
    ax.set_ylim(bottom = 0.0)
    ax.axvline(x = threshold_expected, ls = "--", c = "b", linewidth=1)
    ax.axvline(x = threshold_guide, ls = "--", c = "r", linewidth=1)
    if threshold_auto is not None:
        ax.axvline(x = threshold_auto, ls = "--", c = "g", linewidth=1)
    ax.axvline(x = threshold, ls = "--", c = "k", linewidth=1)
    ax.set_title('KDE of simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Density')

    ax = axes[1, 1]
    ax.plot(sim_x, curv, '-', c='k', lw = 1)
    ax.set_yscale("symlog")
    ax.axvline(x = threshold_expected, ls = "--", c = "b", linewidth=1)
    ax.axvline(x = threshold_guide, ls = "--", c = "r", linewidth=1)
    if threshold_auto is not None:
        ax.axvline(x = threshold_auto, ls = "--", c = "g", linewidth=1)
    ax.axvline(x = threshold, ls = "--", c = "k", linewidth=1)
    ax.set_title('Curvature of simulated doublets')
    ax.set_xlabel('Doublet score')
    ax.set_ylabel('Curvature')

    fig.tight_layout()
    return fig


@timer(logger=logger)
def _run_doublet_detection(
    data: Union[MultimodalData, UnimodalData],
    method: Optional[str] = 'improved',
    raw_mat_key: Optional[str] = 'counts',
    name: Optional[str] = '',
    K: Optional[int] = 100,
    expected_doublet_rate: Optional[float] = None,
    sim_doublet_ratio: Optional[float] = 2.0,
    n_jobs: Optional[int] = -1,
    random_state: Optional[int] = 0,
    plot_hist: Optional[bool] = True,
    manual_correction: Optional[str] = None,
) -> Union[None, Figure]:
    """Calculate doublet scores using a strategy inspired by Scrublet [Wolock18]_for the current data.X; determine a right threshold based on the KDE curve.
       This function should be called after highly_variable_gene selection.

    Parameters
    -----------
    data: ``Union[MultimodalData, UnimodalData]`` object.
        Annotated data matrix with rows for cells and columns for genes. Data must be low quality cell and gene filtered and log-transformed.

    method: ``str``, optional, default: ``improved``
        Choose between ``improved`` and ``scrublet``. If scrublet, assume n_prin_comp=30, k=round(0.5 * sqrt(number of observed cells)), and total neighbors k_adj=round(k * (1.0 + sim_doublet_ratio)).

    raw_mat_key: ``str``, optional, default: ``counts``
        Matrix key for the raw count matrix.

    name: ``str``, optional, default: ``''``
        Name of the sample.

    K: ``int``, optioanl, default: ``100``
        Number of nearest neighbors to query from the observed and simulated cells joint embedding, including the cell itself. pg.neighbors will be called with this parameter. Note K is only used by 'improved' method.

    expected_doublet_rate: ``float``, optional, default: ``None``
        The expected doublet rate for the experiment. By default, calculate the expected rate based on number of cells from the 10x multiplet rate table

    sim_doublet_ratio: ``float``, optional, default: ``2.0``
        The ratio between synthetic doublets and observed cells.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

    random_state: ``int``, optional, default: ``0``
        Random state for doublet simulation, PCA and approximate nearest neighbor search.

    plot_hist: ``bool``, optional, default: ``True``
        If True, plot diagnostic histograms. Each sample would have a figure consisting of 4 panels showing histograms of doublet scores for observed cells (panel 1, density in log scale), simulated doublets (panel 2, density in log scale), KDE plot (panel 3) and signed curvature plot (panel 4, in symlog scale) of doublet scores for simulated doublets.

    manual_correction: ``str``, optional, default: ``None``
        If present, use human guide provided in manual_correction to select threshold. Currently support 'peak', 'expected' and threshold. 'peak' means cutting at the center of the peak and 'expected' means cutting at the expected doublet rate. If not both, convert guide to float and use as user-specified threshold.

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
    >>> pg._run_doublet_detection(data)
    """
    if "highly_variable_features" not in data.var:
        raise ValueError("_run_doublet_detection must be run after highly_variable_features is called!")

    r = sim_doublet_ratio

    if expected_doublet_rate is None:
        expected_doublet_rate = _calc_expected_doublet_rate(data.shape[0])
    rho = expected_doublet_rate

    n_jobs = eff_n_jobs(n_jobs)

    if method == 'improved':
        # Generate simulated doublets data
        datsim = _gen_sim_data(data, raw_mat_key, sim_doublet_ratio, random_state)

        # Concatenate real and simulated data and compute doublet scores
        dat_concat = _concat_and_process(data, datsim, K, r, rho, n_jobs)

        obs_scores = dat_concat.obs.loc[data.obs_names, "doublet_score"].values
        sim_scores = dat_concat.obs.loc[datsim.obs_names, "doublet_score"].values
        obs_pca = data.obsm['X_pca']
    else:
        obs_scores, sim_scores, obs_pca = _run_scrublet(data, raw_mat_key, r, rho, n_jobs, random_state)

    # Calculate theoretical doublet threshold from simulated doublets
    kmeans_labels_, d_emb, d_neo = _calc_expected_emb_rate_in_sim(obs_pca, rho, n_jobs, random_state)
    data.obs["dbl_kmeans_"] = pd.Categorical(kmeans_labels_)
    threshold_guide = np.percentile(sim_scores, d_emb * 100.0 + 1e-6)

    # Calculate expected doublet rate threshold
    threshold_expected = np.percentile(obs_scores, (1.0 - rho) * 100.0 + 1e-6)

    # Determine a scrublet score threshold
    threshold, threshold_auto, x, y, curv = _find_score_threshold(sim_scores, d_neo, threshold_guide, threshold_expected, manual_correction)

    data.obs["doublet_score"] = obs_scores.astype(np.float32)
    data.obs["pred_dbl"] = obs_scores > threshold
    data.uns["doublet_threshold"] = float(threshold)

    neo_dbl_rate = data.obs['pred_dbl'].sum() / data.shape[0]
    neo_sim_dbl_rate = (sim_scores > threshold).sum() / sim_scores.size

    logger.info(f"Sample {name}: total cells = {data.shape[0]}; bimodality coefficient in observed cells = {_calc_bc_sarle(obs_scores):.4f}; expected doublet rate = {rho:.2%}; neotypic doublet rate = {neo_dbl_rate:.2%}; doublet threshold = {threshold:.4f}; bimodality coefficient in simulation = {_calc_bc_sarle(sim_scores):.4f}; neotypic doublet rate in simulation = {neo_sim_dbl_rate:.2%}.")

    fig = None
    if plot_hist:
        fig = _plot_hist(obs_scores, sim_scores, threshold, threshold_guide, threshold_expected, x, y, curv, threshold_auto=threshold_auto)
    return fig


def _identify_doublets_fisher(cluster_labels: Union[pd.Categorical, List[int]], pred_dbl: List[bool], alpha: float = 0.05) -> pd.DataFrame:
    df = pd.crosstab(cluster_labels, pred_dbl)

    if df.shape[1] == 1: # either no doublets or all doublets
        result = pd.DataFrame({'cluster': df.index})
        result['percentage'] = 100.0 if (True in df.columns) else 0.0
        result['pval'] = 1.0
        result['qval'] = 1.0
        return result

    ndbl = df[True].sum().astype(np.int32)
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
    method: Optional[str] = "improved",
    channel_attr: Optional[str] = None,
    clust_attr: Optional[str] = None,
    raw_mat_key: Optional[str] = None,
    min_cell: Optional[int] = 100,
    expected_doublet_rate: Optional[float] = None,
    sim_doublet_ratio: Optional[float] = 2.0,
    K: Optional[int] = 100,
    n_jobs: Optional[int] = -1,
    alpha: Optional[float] = 0.05,
    random_state: Optional[int] = 0,
    plot_hist: Optional[str] = "sample",
    manual_correction: Optional[str] = None,
) -> None:
    """Infer doublets by first calculating doublet scores using a strategy inspired by Scrublet [Wolock18]_ and then smartly determining an appropriate doublet score cutoff [Li20-2]_ .

    This function should be called after clustering if clust_attr is not None. In this case, we will test if each cluster is significantly enriched for doublets using Fisher's exact test.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    method: ``str``, optional, default: ``scrublet``
        Choose between ``scrublet`` and ``improved``.

    channel_attr: ``str``, optional, default: None
        Attribute indicating sample channels. If set, calculate scrublet-like doublet scores per channel.

    clust_attr: ``str``, optional, default: None
        Attribute indicating cluster labels. If set, estimate proportion of doublets in each cluster and statistical significance.

    raw_mat_key: ``str``, optional, default: None
        The key for raw count matrix. By default, Pegasus will first try "counts" and then try "raw.X"

    min_cell: ``int``, optional, default: 100
        Minimum number of cells per sample to calculate doublet scores. For samples having less than 'min_cell' cells, doublet score calculation will be skipped.

    expected_doublet_rate: ``float``, optional, default: ``None``
        The expected doublet rate for the experiment. By default, calculate the expected rate based on number of cells from the 10x multiplet rate table

    sim_doublet_ratio: ``float``, optional, default: ``2.0``
        The ratio between synthetic doublets and observed cells.

    K: ``int``, optional, default: ``100``
        Number of neighbors to query the joint observed & simulated cell embedding.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores.

    alpha: ``float``, optional, default: ``0.05``
        FDR significant level for cluster-level fisher exact test.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    plot_hist: ``str``, optional, default: ``sample``
        If not None, plot diagnostic histograms using ``plot_hist`` as the prefix. If `channel_attr` is None, ``plot_hist.dbl.png`` is generated; Otherwise, ``plot_hist.channel_name.dbl.png`` files are generated. Each figure consists of 4 panels showing histograms of doublet scores for observed cells (panel 1, density in log scale), simulated doublets (panel 2, density in log scale), KDE plot (panel 3) and signed curvature plot (panel 4) of log doublet scores for simulated doublets. Each plot contains two dashed lines. The red dashed line represents the theoretical cutoff (calucalted based on number of cells and 10x doublet table) and the black dashed line represents the cutof inferred from the data.

    manual_correction: ``str``, optional, default: ``None``
        Use human guide to correct doublet threshold for certain channels. This is string representing a comma-separately list. Each item in the list represent one sample and the sample name and correction guide are separated using ':'. The correction guides supported are 'peak', 'expected' and threshold. 'peak' means cutting at the center of the peak; 'expected' means cutting at the expected doublet rate; threshold is the user-specified doublet threshold; if the guide is neither 'peak' nor 'expected', pegasus will try to convert the string into float and use it as doublet threshold. If only one sample available, no need to specify sample name.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs['pred_dbl']``: Predicted singlet/doublet types.

        * ``data.uns['pred_dbl_cluster']``: Only generated if 'clust_attr' is not None. This is a dataframe with two columns, 'Cluster' and 'Qval'. Only clusters with significantly more doublets than expected will be recorded here.

    Examples
    --------
    >>> pg.infer_doublets(data, channel_attr = 'Channel', clust_attr = 'Annotation')
    """
    assert data.get_modality() == "rna"

    if raw_mat_key is None:
        raw_mat_key = 'counts'
        if raw_mat_key not in data.list_keys():
            raw_mat_key = 'raw.X'
    try:
        rawX = data.get_matrix(raw_mat_key)
    except ValueError:
        raise ValueError(f"Cannot detect the raw count matrix {raw_mat_key}; stop inferring doublets!")

    if_plot = plot_hist is not None

    mancor = {}
    if manual_correction is not None:
        if channel_attr is None:
            mancor[''] = manual_correction
        else:
            for item in manual_correction.split(','):
                name, action = item.split(':')
                mancor[name] = action

    if channel_attr is None:
        if data.shape[0] >= min_cell:
            fig = _run_doublet_detection(data, method = method, raw_mat_key = raw_mat_key, expected_doublet_rate = expected_doublet_rate, sim_doublet_ratio = sim_doublet_ratio, \
                                K = K, n_jobs = n_jobs, random_state = random_state, plot_hist = if_plot, manual_correction = mancor.get('', None))
            if if_plot:
                fig.savefig(f"{plot_hist}.dbl.png")
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
                                       {raw_mat_key: rawX[idx]},
                                       {"genome": genome, "modality": modality},
                                       cur_matrix = raw_mat_key)
                # Identify robust genes, count and log normalized and select top 2,000 highly variable features
                identify_robust_genes(unidata)
                log_norm(unidata)
                highly_variable_features(unidata)
                # Run _doublet_detection
                fig = _run_doublet_detection(unidata, method = method, raw_mat_key = raw_mat_key, name = channel, expected_doublet_rate = expected_doublet_rate, sim_doublet_ratio = sim_doublet_ratio, \
                                    K = K, n_jobs = n_jobs, random_state = random_state, plot_hist = if_plot, manual_correction = mancor.get(channel, None))
                if if_plot:
                    fig.savefig(f"{plot_hist}.{channel}.dbl.png")

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





### For old Scrublet implementation, might be obsoleted

def _comp_pca_obs(data, raw_mat_key, pca, n_jobs):
    # subset the raw count matrix
    rawX = data.get_matrix(raw_mat_key)
    obs_umis = rawX.sum(axis = 1, dtype = np.int32).A1
    rawX = rawX[:, data.var["highly_variable_features"].values]

    # standardize and calculate PCA for rawX
    obsX = rawX.astype(np.float32).toarray()
    obsX /= obs_umis.reshape(-1, 1) # normalize each cell

    m1 = obsX.mean(axis = 0) # calculate mean and std
    psum = np.multiply(obsX, obsX).sum(axis=0)
    std = ((psum - obsX.shape[0] * (m1 ** 2)) / (obsX.shape[0] - 1.0)) ** 0.5
    std[std == 0] = 1

    obsX -= m1 # standardize
    obsX /= std

    # compute PCA
    with threadpool_limits(limits = n_jobs):
        obs_pca = pca.fit_transform(obsX.astype(np.float64)) # float64 for reproducibility
        obs_pca = np.ascontiguousarray(obs_pca, dtype=np.float32)

    return rawX, obs_umis, m1, std, obs_pca

def _sim_and_comp_pca(rawX, obs_umis, m1, std, r, random_state, pca, n_jobs):
    # Simulate synthetic doublets
    sim_rawX, pair_idx = _simulate_doublets(rawX, r, random_state)
    sim_umis = obs_umis[pair_idx].sum(axis = 1, dtype = np.int32)

    # standardize and calculate PCA for sim_rawX
    simX = sim_rawX.astype(np.float32).toarray()
    simX /= sim_umis.reshape(-1, 1) # normalize each cell

    simX -= m1 # standardize
    simX /= std

    with threadpool_limits(limits = n_jobs):
        sim_pca = pca.transform(simX) # transform to PC coordinates
        sim_pca = np.ascontiguousarray(sim_pca, dtype=np.float32)

    return sim_pca

def _calc_doublet_score(obs_pca, sim_pca, k, r, rho, n_jobs):
    from pegasus.tools import calculate_nearest_neighbors

    # concatenate observed and simulated data
    pc_coords = np.vstack((obs_pca, sim_pca))
    is_doublet = np.repeat(np.array([0, 1], dtype = np.int32), [obs_pca.shape[0], sim_pca.shape[0]])

    # Calculate k nearest neighbors
    k_adj = int(round(k * (1.0 + r)))
    indices, _, _ = calculate_nearest_neighbors(pc_coords, K=k_adj + 1, n_jobs=n_jobs, exact_k=True)

    # Calculate scrublet-like doublet score
    k_d = is_doublet[indices].sum(axis = 1)
    q = (k_d + 1.0) / (k_adj + 2.0) # Equation 5
    doublet_scores = (q * rho / r) / ((1.0 - rho) - q * (1.0 - rho - rho / r)) # Equation 4
    obs_scores = doublet_scores[0:obs_pca.shape[0]]
    sim_scores = doublet_scores[obs_pca.shape[0]:]

    return obs_scores, sim_scores

def _run_scrublet(
    data: Union[MultimodalData, UnimodalData],
    raw_mat_key: str,
    r: float,
    rho: float,
    n_jobs: int,
    random_state: int,
    n_prin_comps: Optional[int] = 30,
    k: Optional[int] = None,
):
    if k is None:
        k = int(round(0.5 * np.sqrt(data.shape[0])))

    # Compute PC space for original data
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_prin_comps, random_state=random_state)
    rawX, obs_umis, m1, std, obs_pca = _comp_pca_obs(data, raw_mat_key, pca, n_jobs)

    # Simulate synthetic doublets and project to PC space
    sim_pca = _sim_and_comp_pca(rawX, obs_umis, m1, std, r, random_state, pca, n_jobs)

    # Calculatte Scrublet-like doublet scores
    obs_scores, sim_scores = _calc_doublet_score(obs_pca, sim_pca, k, r, rho, n_jobs)

    return obs_scores, sim_scores, obs_pca
