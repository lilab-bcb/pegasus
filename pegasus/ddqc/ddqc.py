import numpy as np
import pandas as pd

from typing import Tuple, Union
from pegasusio import UnimodalData, MultimodalData
from pegasus import tools


# MAD for numpy array. constant is the consistency constant as discussed in https://en.wikipedia.org/wiki/Median_absolute_deviation
# Maybe consider to use scipy.stats.median_abs_deviation(scale = 'normal') to replace?
def _mad(x: np.ndarray, constant: float = 1.4826) -> float:
    return constant * np.median(np.absolute(x - np.median(x)))

# calculates percent ribo for data
def _calculate_percent_ribo(unidata: UnimodalData, ribo_prefix) -> None:
    import re
    ribo_genes = unidata.var_names.map(lambda x: re.match(ribo_prefix, x, flags=re.IGNORECASE)).values.nonzero()[0] # get all genes that match the pattern
    unidata.obs["percent_ribo"] = (unidata.X[:, ribo_genes].sum(axis=1).A1 / np.maximum(unidata.obs["n_counts"].values, 1.0)) * 100 # calculate percent ribo

# function that performs initial qc & clustering; does dimensional reductions and finds DE genes if specified (?)
def _cluster_data(unidata: UnimodalData, n_genes, percent_mito, mito_prefix, ribo_prefix, norm_count: float = 1e5, n_components: int = 50, K: int = 20, resolution: float = 1.3, random_state: int = 29) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    tools.qc_metrics(unidata, min_genes=n_genes, mito_prefix=mito_prefix, percent_mito=percent_mito) # default PG filtering with custom cutoffs
    _calculate_percent_ribo(unidata, ribo_prefix) # calculate percent ribo
    tools.filter_data(unidata) # filtering based on the parameters from qc_metrics
    tools.identify_robust_genes(unidata)

    obs_copy = unidata.obs.copy()
    var_copy = unidata.var.copy()
    uns_copy = unidata.uns.mapping.copy()
    
    tools.log_norm(unidata, norm_count=norm_count)
    tools.highly_variable_features(unidata, consider_batch=False)
    tools.pca(unidata, n_components=n_components, random_state=random_state)
    tools.neighbors(unidata, K=K, random_state=random_state)
    tools.louvain(unidata, resolution=resolution, random_state=random_state)

    return obs_copy, var_copy, uns_copy

# function that performs filtering on the specified metric
# data needs to be clustered
# method - method name for filtering (mad, outlier, cutoff)
# param - parameter for the selected method
# metric name - name of the metric (must be in adata.obs)
# do_upper_co and do_lower_co - whether to do upper and lower cutoff
# record_path - path for recording filtered cells CSVs (keep it None if not needed)
# Only keep "mad" and "outlier"
def _metric_filter(unidata: UnimodalData, method: str, param: float, metric_name: str, do_lower_co: bool = False, do_upper_co: bool = False) -> np.ndarray:
    qc_pass = np.zeros(unidata.shape[0], dtype = bool) # T/F array to tell whether the cell is filtered

    for cl in unidata.obs["louvain_labels"].cat.categories: # iterate though all clusters
        idx = unidata.obs["louvain_labels"] == cl
        values = unidata.obs.loc[idx, metric_name]

        if method == "mad":  # calculate MAD cutoffs, which are median ± param * MAD
            median_v = np.median(values)
            mad_v = _mad(values)
            lower_co = median_v - param * mad_v
            upper_co = median_v + param * mad_v
        else:
            assert method == "outlier" # calculate Outlier cutoffs, which are Q1 - 1.5 * IQR or Q3 + 1.5 * IQR
            q75, q25 = np.percentile(values, [75, 25])
            lower_co = q25 - 1.5 * (q75 - q25)
            upper_co = q75 + 1.5 * (q75 - q25)
        
        qc_pass_cl = np.ones(values.size, dtype = bool)
        if do_lower_co:
            qc_pass_cl &= (values >= lower_co)
        if do_upper_co:
            qc_pass_cl &= (values <= upper_co)
        qc_pass[idx] = qc_pass_cl

    return qc_pass

def _reverse_to_raw_matrix(unidata: UnimodalData, obs_copy: pd.DataFrame, var_copy: pd.DataFrame, uns_copy: dict):
    unidata.obs = obs_copy
    unidata.var = var_copy
    X = unidata.matrices.pop("raw.X")
    unidata.matrices["X"] = X
    unidata.obsm.clear()
    unidata.varm.clear()
    unidata.uns = uns_copy


# function that computes ddqc metrices
# method - method name for filtering (mad or outlier)
# threshold - parameter for the selected method
# do_metric - set to true, if you want to filter the data based on metric
# record_path - path for recording filtered cells CSVs (keep it None if not needed)
def ddqc_metrics(data: Union[MultimodalData, UnimodalData], res=1.3, method="mad", threshold=2, basic_n_genes=100, basic_percent_mito=80, mito_prefix="MT-",
                 ribo_prefix="^Rp[sl]\d", do_counts=True, do_genes=True, do_mito=True, do_ribo=True, random_state=29) -> None:
    if isinstance(data, MultimodalData):
        data = data.current_data()

    obs_copy, var_copy, uns_copy = _cluster_data(data, basic_n_genes, basic_percent_mito, mito_prefix, ribo_prefix, random_state=random_state)

    passed_qc = np.ones(data.shape[0], dtype = bool)
    # for each metric if do_metric is true, the filtering will be performed
    if do_counts:
        passed_qc &= _metric_filter(data, method, threshold, "n_counts", do_lower_co=True)
    if do_genes:
        passed_qc &= _metric_filter(data, method, threshold, "n_genes", do_lower_co=True)
    if do_mito:
        passed_qc &= _metric_filter(data, method, threshold, "percent_mito", do_upper_co=True)
    if do_ribo:
        passed_qc &= _metric_filter(data, method, threshold, "percent_ribo", do_upper_co=True)

    _reverse_to_raw_matrix(data)
    data.obs["passed_qc"] = passed_qc
