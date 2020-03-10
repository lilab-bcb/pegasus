# This code is used to write h5ad 'r+' modified attributes back to disk. The code is adapted from AnnData codebase

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype, is_string_dtype
import h5py
from collections.abc import Mapping
from natsort import natsorted
from scipy.sparse import issparse
import anndata
from typing import List, Tuple



# string dtype
sdt = h5py.string_dtype(encoding = 'utf-8')



def _df_to_records_fixed_width(df: pd.DataFrame) -> Tuple[np.recarray, dict]:
    names = ['index']
    arrays = [df.index.values.astype(sdt)]
    uns_cat = {}
    for col in df.columns:
        names.append(col)
        values = df[col].values
        if is_categorical_dtype(values):
            uns_cat[col + '_categories'] = values.categories
            arrays.append(values.codes)
        elif is_string_dtype(values):
            keywords = np.unique(values)
            if keywords.size < values.size:
                values = pd.Categorical(values, categories = natsorted(keywords))
                uns_cat[col + '_categories'] = values.categories
                arrays.append(values.codes)
            else:
                arrays.append(values.astype(sdt))
        else:
            arrays.append(values)
    formats = [values.dtype for values in arrays]
    return np.rec.fromarrays(arrays, dtype = {'names' : names, 'formats' : formats}), uns_cat


def _to_dict_fixed_width_arrays(data: anndata) -> dict:
    obs_rec, uns_obs = _df_to_records_fixed_width(data.obs)
    var_rec, uns_var = _df_to_records_fixed_width(data.var)
    res = {
        'X': data.X,
        'obs': obs_rec,
        'var': var_rec,
        'obsm': data.obsm,
        'varm': data.varm,
        'uns': {**data.uns, **uns_obs, **uns_var}
    }
    return res



def _parse_whitelist(whitelist: List[str]) -> dict:
    parse_results = {}
    for value in whitelist:
        tokens = value.split("/")
        curr_dict = parse_results
        for i in range(len(tokens) - 1):
            if tokens[i] not in curr_dict:
                curr_dict[tokens[i]] = dict()
            curr_dict = curr_dict[tokens[i]]
            if curr_dict is None:
                break
        if curr_dict is not None:
            curr_dict[tokens[-1]] = None
    return parse_results


def _update_backed_h5ad(group: "hdf5 group", dat: dict, whitelist: dict):
    for key, value in dat.items():
        if not isinstance(key, str):
            logging.warning(
                "Dictionary key {} is transformed to str upon writing to h5,"
                "using string keys is recommended".format(key)
            )
            key = str(key)

        if whitelist is None or key in whitelist:
            if isinstance(value, Mapping):
                subgroup = (
                    group[key] if key in group.keys() else group.create_group(key)
                )
                assert isinstance(subgroup, h5py.Group)
                _update_backed_h5ad(
                    subgroup, value, whitelist[key] if whitelist is not None else None
                )
            else:
                if key in group.keys():
                    del group[key]
                if issparse(value):
                    sparse_mat = group.create_group(key)
                    sparse_mat.attrs["h5sparse_format"] = value.format
                    sparse_mat.attrs["h5sparse_shape"] = np.array(value.shape)
                    sparse_mat.create_dataset(
                        "data", data=value.data, compression="gzip"
                    )
                    sparse_mat.create_dataset(
                        "indices", data=value.indices, compression="gzip"
                    )
                    sparse_mat.create_dataset(
                        "indptr", data=value.indptr, compression="gzip"
                    )
                else:
                    value = np.array(value) if np.ndim(value) > 0 else np.array([value])
                    if value.dtype.kind in {"U", "O"}:
                        value = value.astype(sdt)
                    if value.dtype.names is not None:
                        new_dtype = value.dtype.descr
                        convert_type = False
                        for i in range(len(value.dtype)):
                            if value.dtype[i].kind in {"U", "O"}:
                                new_dtype[i] = (new_dtype[i][0], sdt)
                                convert_type = True
                        if convert_type:
                            value = value.astype(new_dtype)
                    group.create_dataset(key, data=value, compression="gzip")
