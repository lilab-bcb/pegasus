import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from scipy.sparse import issparse, csr_matrix
from typing import Union, List, Tuple, Dict
from anndata import AnnData
from pegasusio import UnimodalData, MultimodalData

import logging
logger = logging.getLogger(__name__)



def eff_n_jobs(n_jobs: int) -> int:
    """ If n_jobs < 0, set it as the number of physical cores _cpu_count """
    if n_jobs > 0:
        return n_jobs

    import psutil
    _cpu_count = psutil.cpu_count(logical=False)
    if _cpu_count is None:
        _cpu_count = psutil.cpu_count(logical=True)
    return _cpu_count


def update_rep(rep: str) -> str:
    """ If rep is None, return rep as mat, which refers to the whole expression matrix
    """
    return rep if rep is not None else "mat"


def X_from_rep(data: AnnData, rep: str, n_comps: int = None) -> np.array:
    """
    If rep is not mat, first check if X_rep is in data.obsm. If not, raise an error.
    If rep is None, return data.X as a numpy array
    """
    if rep != "mat":
        rep_key = "X_" + rep
        if rep_key not in data.obsm.keys():
            raise ValueError("Cannot find {0} matrix. Please run {0} first".format(rep))
        if n_comps is None:
            return data.obsm[rep_key]
        else:
            assert n_comps > 0
            n_dim = min(n_comps, data.obsm[rep_key].shape[1])
            return data.obsm[rep_key][:, 0:n_dim]
    else:
        return data.X if not issparse(data.X) else data.X.toarray()


def W_from_rep(
    data: Union[MultimodalData, UnimodalData, AnnData],
    rep: str,
) -> csr_matrix:
    """
    Return affinity matrix W based on representation rep.
    """
    rep_key = "W_" + rep
    if rep_key not in data.obsp:
        raise ValueError("Affinity matrix does not exist. Please run neighbors first!")
    return data.obsp[rep_key]


# slicing is not designed to work at extracting one element, convert to dense matrix
def slicing(X: Union[csr_matrix, np.ndarray], row: Union[List[bool], List[int], int] = slice(None), col: Union[List[bool], List[int], int] = slice(None), copy: bool = False, squeeze: bool = True) -> np.ndarray:
    result = X[row, col]
    if issparse(X):
        result = result.toarray()
    elif copy:
        result = result.copy()
    if squeeze:
        result = np.squeeze(result)
        if result.ndim == 0:
            result = result.item()
    return result


def calc_mean(X: Union[csr_matrix, np.ndarray], axis: int) -> np.ndarray:
    if not issparse(X):
        return X.mean(axis = axis, dtype = np.float64)

    from pegasus.cylib.fast_utils import calc_mean_sparse
    return calc_mean_sparse(X.shape[0], X.shape[1], X.data, X.indices, X.indptr, axis)


def calc_mean_and_var(X: Union[csr_matrix, np.ndarray], axis: int) -> Tuple[np.ndarray, np.ndarray]:
    if issparse(X):
        from pegasus.cylib.fast_utils import calc_mean_and_var_sparse
        return calc_mean_and_var_sparse(X.shape[0], X.shape[1], X.data, X.indices, X.indptr, axis)
    else:
        from pegasus.cylib.fast_utils import calc_mean_and_var_dense
        return calc_mean_and_var_dense(X.shape[0], X.shape[1], X, axis)


def calc_expm1(X: Union[csr_matrix, np.ndarray]) -> np.ndarray:
    if not issparse(X):
        return np.expm1(X)
    res = X.copy()
    np.expm1(res.data, out = res.data)
    return res


def calc_stat_per_batch(X: Union[csr_matrix, np.ndarray], batch: Union[pd.Categorical, np.ndarray, list]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    from pandas.api.types import is_categorical_dtype
    if is_categorical_dtype(batch):
        nbatch = batch.categories.size
        codes = batch.codes.astype(np.int32)
    else:
        codes = np.array(batch, dtype = np.int32)
        nbatch = codes.max() + 1 # assume cluster label starts from 0

    if issparse(X):
        from pegasus.cylib.fast_utils import calc_stat_per_batch_sparse
        return calc_stat_per_batch_sparse(X.shape[0], X.shape[1], X.data, X.indices, X.indptr, nbatch, codes)
    else:
        from pegasus.cylib.fast_utils import calc_stat_per_batch_dense
        return calc_stat_per_batch_dense(X.shape[0], X.shape[1], X, nbatch, codes)


def normalize_by_count(X: Union[csr_matrix, np.ndarray], robust: List[bool], norm_count: float, log_transform: bool) -> np.ndarray:
    scale = None
    if issparse(X):
        from pegasus.cylib.fast_utils import normalize_by_count_sparse
        scale = normalize_by_count_sparse(X.shape[0], X.shape[1], X.data, X.indices, X.indptr, robust, norm_count)
        if log_transform:
            np.log1p(X.data, out = X.data)
    else:
        from pegasus.cylib.fast_utils import normalize_by_count_dense
        scale = normalize_by_count_dense(X.shape[0], X.shape[1], X, robust, norm_count)
        if log_transform:
            np.log1p(X, out = X)
    return scale


def calc_sig_background(X: Union[csr_matrix, np.ndarray], bins: pd.Categorical, mean_vec: List[float]) -> Tuple[np.ndarray, np.ndarray]:
    n_bins = bins.categories.size
    codes = bins.codes.astype(np.int32)

    if mean_vec.dtype == np.float32:
        mean_vec = mean_vec.astype(np.float64)
    if X.dtype == np.float64: # not sure if we should have this; we may consider make all matrices float32?
        X = X.astype(np.float32)

    if issparse(X):
        from pegasus.cylib.fast_utils import calc_sig_background_sparse
        return calc_sig_background_sparse(X.shape[0], X.shape[1], X.data, X.indices, X.indptr, n_bins, codes, mean_vec)
    else:
        from pegasus.cylib.fast_utils import calc_sig_background_dense
        return calc_sig_background_dense(X.shape[0], X.shape[1], X, n_bins, codes, mean_vec)


def simulate_doublets(X: Union[csr_matrix, np.ndarray], sim_doublet_ratio: float, random_state: int = 0) -> Tuple[Union[csr_matrix, np.ndarray], np.ndarray]:
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


def check_batch_key(data: Union[MultimodalData, UnimodalData], batch: str, warning_msg: str) -> bool:
    if batch is None:
        return False

    if batch not in data.obs:
        logger.warning(f"Batch key {batch} does not exist. {warning_msg}")
        return False
    else:
        if not is_categorical_dtype(data.obs[batch]):
            data.obs[batch] = pd.Categorical(data.obs[batch].values)
        if data.obs[batch].cat.categories.size == 1:
            logger.warning(f"Batch key {batch} only contains one batch. {warning_msg}")
            return False
    return True



import pkg_resources

predefined_signatures = dict(
    cell_cycle_human=pkg_resources.resource_filename("pegasus", "data_files/cell_cycle_human.gmt"),
    cell_cycle_mouse=pkg_resources.resource_filename("pegasus", "data_files/cell_cycle_mouse.gmt"),
    gender_human=pkg_resources.resource_filename("pegasus", "data_files/gender_human.gmt"),
    gender_mouse=pkg_resources.resource_filename("pegasus", "data_files/gender_mouse.gmt"),
    mitochondrial_genes_human=pkg_resources.resource_filename("pegasus", "data_files/mitochondrial_genes_human.gmt"),
    mitochondrial_genes_mouse=pkg_resources.resource_filename("pegasus", "data_files/mitochondrial_genes_mouse.gmt"),
    ribosomal_genes_human=pkg_resources.resource_filename("pegasus", "data_files/ribosomal_genes_human.gmt"),
    ribosomal_genes_mouse=pkg_resources.resource_filename("pegasus", "data_files/ribosomal_genes_mouse.gmt"),
    apoptosis_human=pkg_resources.resource_filename("pegasus", "data_files/apoptosis_human.gmt"),
    apoptosis_mouse=pkg_resources.resource_filename("pegasus", "data_files/apoptosis_mouse.gmt"),
)

predefined_pathways = dict(
    hallmark=pkg_resources.resource_filename("pegasus", "data_files/h.all.v7.5.1.symbols.gmt"),
    canonical_pathways=pkg_resources.resource_filename("pegasus", "data_files/c2.cp.v7.5.1.symbols.gmt"),
)

def load_signatures_from_file(input_file: str) -> Dict[str, List[str]]:
    signatures = {}
    with open(input_file) as fin:
        for line in fin:
            items = line.strip().split('\t')
            signatures[items[0]] = list(set(items[2:]))
    logger.info(f"Loaded signatures from GMT file {input_file}.")
    return signatures


def largest_variance_from_random_matrix(
    ncells: int,
    nfeatures: int,
    pval: str = "0.05",
) -> float:
    """ Select the largest variance from a random generated matrix. See [Johnstone 2001](https://projecteuclid.org/journals/annals-of-statistics/volume-29/issue-2/On-the-distribution-of-the-largest-eigenvalue-in-principal/10.1214/aos/1009210544.full) and [Shekhar et al. 2022](https://elifesciences.org/articles/73809) for more details.

    Parameters
    ----------
    ncells: ``int``.
        Number of cells.

    nfeatures: ``int``.
        Number of featuers (e.g. highly variable genes)

    pval: ``str``, optional (default: "0.05").
        P value cutoff on the null distribution (random matrix), choosing from "0.01" and "0.05".

    Returns
    -------
    ``float``, the largest variance from a random matrix at pval.

    Examples
    --------
    >>> pg.largest_variance_from_random_matrix(10000, 2000)
    """
    quantiles = {"0.01": 2.023335, "0.05": 0.9792895} # quantiles from the Tracy-Widom distribution of order 1.
    assert pval in ["0.01", "0.05"]
    val1 = (ncells - 1) ** 0.5
    val2 = nfeatures ** 0.5
    mu = (val1 + val2) ** 2
    sigma = (val1 + val2) * (1.0 / val1 + 1.0 / val2) ** (1.0 / 3.0)
    res = (quantiles[pval] * sigma + mu) / (ncells - 1)

    return res
