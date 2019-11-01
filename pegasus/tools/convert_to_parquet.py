import gc
import json
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse
from pegasus.io import read_input

from .. import decorators as pg_deco

obsm_whitelist = ['X_pca', 'X_rpca', 'X_tsne', 'X_fitsne', 'X_umap', 'X_fle', 'X_net_tsne', 'X_net_umap', 'X_net_fle']
obsm_whitelist_3d = ['X_diffmap_pca']


def create_schema(data):
    df = data.obs.iloc[[0, 1]].copy()
    obs_coords = []
    columns_to_add = []
    obsm = []
    for key in data.obsm.keys():
        ndim = None
        if key in obsm_whitelist:
            ndim = 2
        elif key in obsm_whitelist_3d:
            ndim = 3
        if ndim is not None:
            obs_coords.append(key)
            coordinate_columns = ['{}_{}'.format(key, i) for i in range(1, ndim + 1)]
            columns_to_add += coordinate_columns
            obsm.append({'name': key, 'dimensions': ndim})

    columns_to_add += data.var_names.to_list()
    empty_data = np.zeros((2, len(columns_to_add)), dtype='float32')
    empty_df = pd.DataFrame(index=df.index, columns=columns_to_add, data=empty_data)
    df = df.join(empty_df)
    table = pa.Table.from_pandas(df)
    schema = table.schema

    schema = schema.with_metadata(
        {b'pegasus': json.dumps(
            {'obsm': obsm, 'var': data.var.index.values.tolist(), 'obs': data.obs.columns.values.tolist()}).encode(
            'utf8')})
    return schema


def to_df(data):
    X = data.X
    if scipy.sparse.issparse(X):
        X = X.toarray()

    if str(X.dtype) != 'float32':
        X = X.astype('float32')

    df = pd.DataFrame(index=data.obs.index, data=X, columns=data.var_names)
    for key in data.obsm.keys():
        if key in obsm_whitelist:
            df["{}_1".format(key)] = data.obsm[key][:, 0].astype('float32')
            df["{}_2".format(key)] = data.obsm[key][:, 1].astype('float32')
        elif key in obsm_whitelist_3d:
            df["{}_1".format(key)] = data.obsm[key][:, 0].astype('float32')
            df["{}_2".format(key)] = data.obsm[key][:, 1].astype('float32')
            df["{}_3".format(key)] = data.obsm[key][:, 2].astype('float32')
    df = df.join(data.obs)
    return df

@pg_deco.TimeLogger()
def convert_to_parquet(data, output_name, nthreads, row_group_size):
    if not output_name.endswith(".pq") and not output_name.endswith(".parquet"):
        output_name = output_name + '.pq'
    schema = create_schema(data)
    with pq.ParquetWriter(output_name, schema) as writer:
        for i in range(0, data.shape[0], row_group_size):
            end = i + row_group_size
            end = min(end, data.shape[0])
            df = to_df(data[i:end])
            table = pa.Table.from_pandas(df, schema=schema, nthreads=nthreads)
            writer.write_table(table)
            gc.collect()
    print(output_name + " is written!")


def run_conversion(input_h5ad_file, output_name, nthreads, row_group_size):
    data = read_input(input_h5ad_file)
    convert_to_parquet(data, output_name, nthreads, row_group_size)