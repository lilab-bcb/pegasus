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
def _calculate_percent_ribo(data: MultimodalData, ribo_prefix) -> None:
    import re
    ribo_genes = data.var_names.map(lambda x: re.match(ribo_prefix, x, flags=re.IGNORECASE) is not None).values.nonzero()[0] # get all genes that match the pattern
    data.obs["percent_ribo"] = (data.X[:, ribo_genes].sum(axis=1).A1 / np.maximum(data.obs["n_counts"].values, 1.0)) * 100 # calculate percent ribo

# function that performs initial qc & clustering; does dimensional reductions and finds DE genes if specified (?)
def _cluster_data(data: MultimodalData, n_genes: int, percent_mito: float, mito_prefix: str, ribo_prefix: str, norm_count: float = 1e5, n_components: int = 50, K: int = 20, resolution: float = 1.3, random_state: int = 29) -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    tools.qc_metrics(data, min_genes=n_genes, mito_prefix=mito_prefix, percent_mito=percent_mito) # default PG filtering with custom cutoffs
    _calculate_percent_ribo(data, ribo_prefix) # calculate percent ribo
    tools.filter_data(data) # filtering based on the parameters from qc_metrics
    tools.identify_robust_genes(data)

    obs_copy = data.obs.copy()
    var_copy = data.var.copy()
    var_copy.drop(columns=["robust", "highly_variable_features", "n_cells", "percent_cells"], inplace=True, errors="ignore")
    uns_copy = data.uns.mapping.copy()

    tools.log_norm(data, norm_count=norm_count)
    tools.highly_variable_features(data, consider_batch=False)
    tools.pca(data, n_components=n_components, random_state=random_state)
    tools.neighbors(data, K=K, random_state=random_state)
    tools.louvain(data, resolution=resolution, random_state=random_state)

    return obs_copy, var_copy, uns_copy

# function that performs filtering on the specified metric
# data needs to be clustered
# method - method name for filtering (mad, outlier, cutoff)
# param - parameter for the selected method
# metric name - name of the metric (must be in adata.obs)
# do_upper_co and do_lower_co - whether to do upper and lower cutoff
# record_path - path for recording filtered cells CSVs (keep it None if not needed)
# Only keep "mad" and "outlier"
def _metric_filter(data: MultimodalData, method: str, param: float, metric_name: str, do_lower_co: bool = False, do_upper_co: bool = False) -> np.ndarray:
    qc_pass = np.zeros(data.shape[0], dtype = bool) # T/F array to tell whether the cell is filtered

    for cl in data.obs["louvain_labels"].cat.categories: # iterate though all clusters
        idx = data.obs["louvain_labels"] == cl
        values = data.obs.loc[idx, metric_name]

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
def ddqc_metrics(data: MultimodalData, res=1.3, method="mad", threshold=2, basic_n_genes=100, basic_percent_mito=80, mito_prefix="MT-",
                 ribo_prefix="^Rp[sl]\d", do_counts=True, do_genes=True, do_mito=True, do_ribo=True, random_state=29) -> None:
    assert isinstance(data, MultimodalData)
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

    _reverse_to_raw_matrix(data.current_data(), obs_copy, var_copy, uns_copy)
    data.obs["passed_qc"] = passed_qc
