import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from typing import List, Union
from pegasusio import UnimodalData, MultimodalData

import anndata

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



def _gen_anndata(
    data: Union[MultimodalData, UnimodalData],
    obs_columns: List[str],
    features: str = "highly_variable_features",
    matkey: str = "raw.X",
) -> anndata.AnnData:
    """ Generate a new Anndata object for scvitools

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    obs_columns: ``List[str]``
        A list of obs keys that should be included in the new anndata.
    features: ``str``, optional, default: ``highly_variable_features``
        a keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.
    matkey: ``str``, optional, default: ``raw.X``
        Matrix key for the raw count

    Returns
    -------
    An AnnData object.
    """
    mat = data.get_matrix(matkey)
    obs_field = data.obs[obs_columns]
    if features is not None:
        assert features in data.var
        idx = data.var[features].values
        X = mat[:, idx]
        var_field = data.var.loc[idx, []]
    else:
        X = mat
        var_field = data.var[[]]

    return anndata.AnnData(X = X, obs = obs_field, var = var_field)



@njit(fastmath=True, cache=True)
def _select_csr(data, indices, indptr, indexer, new_size):
    data_new = np.zeros_like(data[0:new_size])
    indices_new = np.zeros_like(indices[0:new_size])
    indptr_new = np.zeros_like(indptr)

    cnt = 0
    for i in range(indptr.size - 1):
        indptr_new[i] = cnt
        for j in range(indptr[i], indptr[i + 1]):
            new_idx = indexer[indices[j]]
            if new_idx >= 0:
                data_new[cnt] = data[j]
                indices_new[cnt] = new_idx
                cnt += 1
    indptr_new[indptr.size - 1] = cnt    

    return data_new, indices_new, indptr_new


def _gen_query_anndata(
    data: Union[MultimodalData, UnimodalData],
    obs_columns: List[str],
    ref_features: pd.Index,
    matkey: str = "raw.X",
) -> anndata.AnnData:
    """ Generate a new query Anndata object for scvitools

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    obs_columns: ``List[str]``
        A list of obs keys that should be included in the new anndata.
    ref_features: ``pd.Index``
        A pandas index of reference feature names
    matkey: ``str``, optional, default: ``raw.X``
        Matrix key for the raw count

    Returns
    -------
    An AnnData object.
    """
    mat = data.get_matrix(matkey)
    obs_field = data.obs[obs_columns]
    var_field = pd.DataFrame(index = ref_features)

    indexer = ref_features.get_indexer(data.var_names)
    new_size = (indexer[mat.indices]>=0).sum()
    data_new, indices_new, indptr_new = _select_csr(mat.data, mat.indices, mat.indptr, indexer, new_size)
    X = csr_matrix((data_new, indices_new, indptr_new), shape = (mat.shape[0], ref_features.size))
    X.sort_indices()

    return anndata.AnnData(X = X, obs = obs_field, var = var_field)
