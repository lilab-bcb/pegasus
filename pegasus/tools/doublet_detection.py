import time
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from sklearn.decomposition import TruncatedSVD
from sklearn.mixture import GaussianMixture

from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection as fdr

from typing import List, Optional, Union

import logging
logger = logging.getLogger(__name__)

from pegasusio import MultimodalData
from pegasusio import timer

from .clustering import partition_cells_by_kmeans


def _one_tail_test(scores: List[float], mean: float, std: float, alpha: float = 0.05) -> List[bool]:
    idx = scores > mean
    pvals = 1.0 - norm.cdf(scores[idx], loc = mean, scale = std)
    passed, qvals = fdr(pvals, alpha = alpha)

    outliers = np.zeros(scores.size, dtype = bool)
    outliers[idx] = passed

    return outliers


def _identify_cell_doublets(scores: List[float], alpha: float = 0.05, random_state: int = 0):
    scores = np.log(scores) # log transformed
    scores_reshaped = scores.reshape(-1, 1)
    # First fit three normal distributions
    gm = GaussianMixture(n_components = 3, random_state = random_state)
    gm.fit(scores_reshaped)
    means = gm.means_.ravel()
    stds = np.sqrt(gm.covariances_.ravel())
    pos = np.argsort(means)[1]
    outliers = _one_tail_test(scores, means[pos], stds[pos], alpha)
    prev_ndbl = outliers.sum()
    # Fit two normals by excluding outliers
    gm.set_params(n_components = 2)
    gm.fit(scores_reshaped[~outliers])
    means = gm.means_.ravel()
    stds = np.sqrt(gm.covariances_.ravel())
    pos = np.argsort(means)[1]
    outliers = _one_tail_test(scores, means[pos], stds[pos], alpha)
    ndbl = outliers.sum()
    # Iteratively reduce false cell doublets
    gm.set_params(warm_start = True)
    while ndbl < prev_ndbl:
        gm.fit(scores_reshaped[~outliers])
        means = gm.means_.ravel()
        stds = np.sqrt(gm.covariances_.ravel())
        pos = np.argsort(means)[1]
        outliers = _one_tail_test(scores, means[pos], stds[pos], alpha)
        prev_ndbl = ndbl
        ndbl = outliers.sum()
    # Predict singlets and transition cells
    preds = gm.predict(scores_reshaped[~outliers])
    if pos == 0:
        preds = 1 - preds
    # Generate labels: 0, singlets; 1, transition; 2, doublets
    labels = np.zeros(scores.size)
    labels[outliers] = 2
    labels[~outliers] = preds

    return labels


def _identify_doublets_fisher(cluster_labels: Union[pd.Categorical, List[int]], dbl_codes: List[int], alpha: float = 0.05) -> pd.DataFrame:
    dbls = dbl_codes > 1
    df = pd.crosstab(cluster_labels, dbls)

    ndbl = df[True].sum()
    a = df[True].values.astype(np.uint)
    b = df[False].values.astype(np.uint)
    c = ndbl - a
    d = (dbl_codes.size - ndbl) - b

    avg_dblr = ndbl / dbl_codes.size
    freqs = a / (a + b)

    import fisher
    pvals = fisher.pvalue_npy(a, b, c, d)[2]
    passed, qvals = fdr(pvals, alpha = alpha)

    posvec = np.where(passed)[0][freqs[passed] > avg_dblr]

    result = pd.DataFrame({'cluster': df.index[posvec], 'percentage': freqs[posvec] * 100.0, 'pval': pvals[posvec], 'qval': qvals[posvec]})
    result.sort_values('percentage', ascending = False, inplace = True)
    result.reset_index(drop=True, inplace=True)

    return result



@timer(logger=logger)
def infer_doublets(
    data: MultimodalData,
    dbl_attr: Optional[str] = 'scrublet_scores',
    channel_attr: Optional[str] = None,
    clust_attr: Optional[str] = None,
    n_components: Optional[int] = 50,
    robust: Optional[bool] = True,
    n_clusters: Optional[int] = 30,
    n_clusters2: Optional[int] = 50,
    n_init: Optional[int] = 10,
    min_avg_cells_per_final_cluster: Optional[int] = 10,
    alpha: Optional[float] = 0.05,
    random_state: Optional[int] = 0,
    verbose: Optional[bool] = False,
) -> None:
    """Infer doublets based on i.e. Scrublet scores.
    Inspired by Pijuan-Sala et al. Nature 2019 and Popescu et al. Nature 2019.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    dbl_attr: ``str``, optional, default: ``scrublet_scores``
        Attribute indicating calculated doublet scores from Scrublet.

    channel_attr: ``str``, optional, default: None
        Attribute indicating sample channels. If None, assuming cell-level and sample-level doublets are already calculated and saved in 'data.obs[pred_dbl_type]'.

    clust_attr: ``str``, optional, default: None
        Attribute indicating cluster labels. If None, does not perform cluster-level doublet detection.

    n_components: ``int``, optional, default: ``50``
        Number of PC components for sample-level doublet inference. Note that we use all genes (sparse matrix) and truncated SVD to infer PCs. Because truncated SVD does not reduce means and the first component correlates with the mean, we will use n_components + 1 in sklearn.decomposition.TruncatedSVD .

    robust: ``bool``, optional, default ``True``
        If robust == True, use algorithm = 'arpack'; otherwise, use algorithm = 'randomized'.

    n_clusters: ``int``, optional, default: ``30``
        The number of first level clusters.

    n_clusters2: ``int``, optional, default: ``50``
        The number of second level clusters.

    n_init: ``int``, optional, default: ``10``
        Number of kmeans tries for the first level clustering. Default is set to be the same as scikit-learn Kmeans function.

    min_avg_cells_per_final_cluster: ``int``, optional, default: ``10``
        Expected number of cells per cluster after the two-level KMeans.

    alpha: ``float``, optional, default: ``0.05``
        FDR significant level for statistical tests.

    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    verbose: ``bool``, optional, default: ``False``
        If true, pegasus generates density plots for each channel under working directory with name channel.dbl.png. In addition, diagnostic outputs will be generated.

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

    if channel_attr is not None:
        assert is_categorical_dtype(data.obs[channel_attr])
        from pegasus.plotting import doublet_plot

        dbl_codes = np.zeros(data.shape[0], dtype = np.int32)

        channels = data.obs[channel_attr].cat.categories
        for channel in channels:
            idx = (data.obs[channel_attr] == channel).values
            dbl_scores = data.obs.loc[idx, dbl_attr].values

            # Cell-level doublet
            dblc_codes = _identify_cell_doublets(dbl_scores, alpha = alpha, random_state = random_state) # dblc doublet at cell level

            idx_dblc = dblc_codes == 2
            idx_dblnc = ~idx_dblc
            ncdbl = idx_dblc.sum()
            freqs = []

            # Sample-level test
            # Truncated SVD including all genes
            if ncdbl > 0:
                tsvd = TruncatedSVD(n_components = n_components + 1, algorithm = 'arpack' if robust else 'randomized', random_state = random_state)
                X_tpca = np.ascontiguousarray(tsvd.fit_transform(data.X[idx]))

                clusters = partition_cells_by_kmeans(X_tpca, n_clusters, n_clusters2, n_init, random_state, min_avg_cells_per_final_cluster)

                sigs = _identify_doublets_fisher(clusters, dblc_codes, alpha = alpha) # significant clusters

                for cluster in sigs['cluster']:
                    idxc = clusters == cluster
                    idx_dbls = idxc & idx_dblnc
                    dblc_codes[idx_dbls] = 3
                    freqs.append(1.0 - idx_dbls.sum() / idxc.sum())

            # assign channel predictions to dbl_codes
            dbl_codes[idx] = dblc_codes

            # QC statistics
            nsdbl = (dblc_codes == 3).sum()
            min_score = dbl_scores[idx_dblc].min() if ncdbl > 0 else None
            min_freq = min(freqs) if len(freqs) > 0 else None

            if verbose:
                fig = doublet_plot(dbl_scores, dblc_codes, return_fig = True)
                fig.savefig(f"{channel}.dbl.png")
                logger.info(f"Channel {channel}: {ncdbl} cell-level doublets and {nsdbl} sample-level doublets were detected!")
                if min_score is not None:
                    logger.info(f"Doublet score cutoff for cell-level duoblets is {min_score:.3f}.")
                if min_freq is not None:
                    logger.info(f"Doublet frequency cutoff for sample-level doublets is {min_freq:.3f}.")
                logger.info(f"Density plot {channel}.dbl.png is generated.")

            logger.info(f"Channel {channel} contains {dbl_scores.size} cells and {ncdbl + nsdbl} predicted doublets. The predicted doublet rate is {(ncdbl+nsdbl)/dbl_scores.size:.2%}.")

        data.obs['pred_dbl_type'] = pd.Categorical.from_codes(dbl_codes, categories=['singlet', 'singlet-2', 'doublet-cell', 'doublet-sample'])

    if clust_attr is not None:
        clusters = data.obs[clust_attr].values
        data.uns['pred_dbl_cluster'] = _identify_doublets_fisher(clusters, dbl_codes, alpha = alpha)

    logger.info('Doublets are predicted!')



def mark_singlets(
    data: MultimodalData,
    demux_attr: Optional[str] = 'demux_type',
    dbl_clusts: Optional[str] = None,
) -> None:
    """Convert doublet prediction into doublet annotations that Pegasus can recognize. Must run infer_doublets first.

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
    >>> pg.mark_singlets(data, dbl_clusts='Annotation:B/T doublets')
    """
    data.obs[demux_attr] = 'singlet'
    idx = data.obs['pred_dbl_type'].map(lambda x: x.startswith('doublet'))
    data.obs.loc[idx, demux_attr] = 'doublet'
    if dbl_clusts is not None:
        cluster, value_str = dbl_clusts.split(':')
        idx = np.isin(data.obs[cluster], value_str.split(','))
        data.obs.loc[idx, demux_attr] = 'doublet'

        codes = data.obs['pred_dbl_type'].values.codes
        idx_4 = idx & (codes < 2)
        if idx_4.sum() > 0:
            codes = codes.copy()
            codes[idx_4] = 4
            data.obs['pred_dbl_type'] = pd.Categorical.from_codes(codes, categories = ['singlet', 'singlet-2', 'doublet-cell', 'doublet-sample', 'doublet-cluster'])
    
    data.obs[demux_attr] = pd.Categorical(data.obs[demux_attr], categories = ['singlet', 'doublet'])
