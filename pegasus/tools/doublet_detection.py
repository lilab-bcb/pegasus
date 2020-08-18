import time
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from sklearn.decomposition import TruncatedSVD

from typing import List, Optional

import logging
logger = logging.getLogger(__name__)

from pegasusio import MultimodalData
from pegasusio import timer

from .clustering import partition_cells_by_kmeans


def _one_tail_test(scores: List[float], alpha: float) -> List[bool]:
    # Calculate median, MAD and std
    median = np.median(scores)
    mad = np.median(scores[scores > median] - median)
    std = 1.4826 * mad # sigma = k dot MAD, k = 1.4826
    # Conduct one-side test because the alternative hypothesis should not have a smaller doublet score.
    from scipy.stats import norm
    from statsmodels.stats.multitest import fdrcorrection as fdr
    pvals = 1.0 - norm.cdf(scores, loc = median, scale = std)
    passed, qvals = fdr(pvals, alpha = alpha)

    return passed


@timer(logger=logger)
def infer_doublets(
    data: MultimodalData,
    clust_attr: str,
    channel_attr: Optional[str] = 'Channel',
    dbl_attr: Optional[str] = 'scrublet_scores',
    n_components: Optional[int] = 50,
    robust: Optional[bool] = True,
    n_clusters: Optional[int] = 30,
    n_clusters2: Optional[int] = 50,
    n_init: Optional[int] = 10,
    min_avg_cells_per_final_cluster: Optional[int] = 20,
    alpha_cell: Optional[float] = 0.05,
    alpha_sample: Optional[float] = 0.05,
    alpha_all: Optional[float] = 0.05,
    random_state: Optional[int] = 0,
) -> None:
    """Infer doublets based on i.e. Scrublet scores.
    Following Pijuan-Sala et al. Nature 2019 and Popescu et al. Nature 2019.
    
    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    clust_attr: ``str``
        Attribute indicating cluster labels.

    channel_attr: ``str``, optional, default: ``Channel``
        Attribute indicating sample channels.

    dbl_attr: ``str``, optional, default: ``scrublet_scores``
        Attribute indicating calculated doublet scores from Scrublet.

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

    alpha_cell: ``float``, optional, default: ``0.05``
        FDR significant level for per cell test.

    alpha_sample: ``float``, optional, default: ``0.05``
        FDR significant level for per sample test.

    alpha_all: ``float``, optional, default: ``0.05``
        FDR significant level for per overall cluster test.
        
    random_state: ``int``, optional, default: ``0``
        Random seed for reproducing results.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs[class_label]``: Cluster labels for cells as categorical data.

    Examples
    --------
    >>> pg.infer_doublets(data)
    """
    assert is_categorical_dtype(data.obs[channel_attr])

    dbl_codes = np.zeros(data.shape[0], dtype = np.int32)

    channels = data.obs[channel_attr].cat.categories
    for channel in channels:
        start = time.perf_counter()
        idx = (data.obs[channel_attr] == channel).values
        dbl_scores = data.obs.loc[idx, dbl_attr].values
        # Cell-level test
        passed_cell = _one_tail_test(dbl_scores, alpha_cell)
        # Sample-level test
        # Truncated SVD including all genes
        tsvd = TruncatedSVD(n_components = n_components + 1, algorithm = 'arpack' if robust else 'randomized', random_state = random_state)
        X_tpca = np.ascontiguousarray(tsvd.fit_transform(data.X[idx]))

        labels = partition_cells_by_kmeans(X_tpca, n_clusters, n_clusters2, n_init, random_state, min_avg_cells_per_final_cluster)
        
        cats = np.unique(labels)
        mdscores = np.zeros(cats.size, dtype = np.float32)
        for i in range(cats.size):
            mdscores[i] = np.median(dbl_scores[labels == cats[i]])

        passed = _one_tail_test(mdscores, alpha_sample)
        passed_sample = np.isin(labels, cats[passed])

        posvec = np.where(idx)[0]
        dbl_codes[posvec[passed_sample]] = 2
        dbl_codes[posvec[passed_cell]] = 1

        end = time.perf_counter()
        logger.info(f"Channel {channel} is processed, time spent = {end-start:.2f}s.")

    cats = data.obs[clust_attr].cat.categories
    freqs = np.zeros(cats.size, dtype = np.float32)
    for i in range(cats.size):
        idx = (data.obs[clust_attr] == cats[i]).values
        freqs[i] = (dbl_codes[idx] > 0).sum() * 1.0 / idx.sum()

    passed = _one_tail_test(freqs, alpha_all)
    passed_all = np.isin(data.obs[clust_attr].values, cats[passed])
    dbl_codes[passed_all & (dbl_codes == 0)] = 3

    data.obs['pred_dbl_type'] = pd.Categorical.from_codes(dbl_codes, categories = ['singlet', 'doublet-cell', 'doublet-sample', 'doublet-cluster'])


