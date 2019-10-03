import numpy as np
from scipy.sparse import issparse


def update_rep(rep: str) -> str:
    """ If rep is None, return rep as mat, which refers to the whole expression matrix
    """
    return rep if rep is not None else "mat"


def X_from_rep(data: "AnnData", rep: str) -> np.array:
    """
    If rep is not mat, first check if X_rep is in data.obsm. If not, raise an error. 
    If rep is None, return data.X as a numpy array
    """
    if rep != "mat":
        rep_key = "X_" + rep
        if rep_key not in data.obsm.keys():
            raise ValueError("Cannot find {0} matrix. Please run {0} first".format(rep))
        return data.obsm[rep_key]
    else:
        return data.X if not issparse(data.X) else data.X.toarray()


def W_from_rep(data: "AnnData", rep: str) -> "csr_matrix":
    """
    Return affinity matrix W based on representation rep.
    """
    rep_key = "W_" + rep
    if rep_key not in data.uns:
        raise ValueError("Affinity matrix does not exist. Please run neighbors first!")
    return data.uns[rep_key]


def knn_is_cached(
    data: "AnnData", indices_key: str, distances_key: str, K: int
) -> bool:
    return (
        (indices_key in data.uns)
        and (distances_key in data.uns)
        and data.uns[indices_key].shape[0] == data.shape[0]
        and (K <= data.uns[indices_key].shape[1] + 1)
    )
