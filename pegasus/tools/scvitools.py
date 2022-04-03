import numpy as np
import pandas as pd

from typing import List, Union
from pegasusio import UnimodalData, MultimodalData

import anndata

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



def _gen_anndata(
    data: Union[MultimodalData, UnimodalData],
    obs_columns: List[str],
    matkey: str = "raw.X",
    features: str = "highly_variable_features",
) -> anndata.AnnData:
    """ Generate a new Anndata object for scvitools

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    obs_columns: ``List[str]``
        A list of obs keys that should be included in the new anndata.
    matkey: ``str``, optional, default: ``raw.X``
        Matrix key for the raw count
    features: ``str``, optional, default: ``highly_variable_features``
        a keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.

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
