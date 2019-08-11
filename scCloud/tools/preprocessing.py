import time
import numpy as np
import pandas as pd
import xlsxwriter

from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.utils.extmath import randomized_svd

from typing import List, Tuple


def qc_metrics(data: 'AnnData', mito_prefix: str ='MT-', min_genes: int = 500,
               max_genes: int = 6000, min_umis: int = 100, max_umis: int = 600000, 
               percent_mito: float = 0.1,
               percent_cells: float = 0.0005) -> None:
    """
    TODO: documentation

    Sets passed_qc, n_genes, n_counts, percent_mito on data.obs and passed_qc, n_cells, percent_cells, and robust on data.var

    :param data:
       Annotated data matrix
    :param mito_prefix: str
       String that mitochrondrial genes start with
    :param percent_cells: float
       Cutoff for a feature to be `robust`
    """

    data.obs['passed_qc'] = False

    data.obs['n_genes'] = data.X.getnnz(axis=1)
    data.obs['n_counts'] = data.X.sum(axis=1).A1

    mito_prefixes = mito_prefix.split(',')

    def startswith(name):
        for prefix in mito_prefixes:
            if name.startswith(prefix):
                return True
        return False

    mito_genes = data.var_names.map(startswith).values.nonzero()[0]
    data.obs['percent_mito'] = data.X[:, mito_genes].sum(axis=1).A1 / np.maximum(
        data.obs['n_counts'].values, 1.0
    )

    # Assign passed_qc
    filters = [data.obs['n_genes'] >= min_genes, data.obs['n_genes'] < max_genes, data.obs['n_counts'] >= min_umis,
               data.obs['n_counts'] < max_umis, data.obs['percent_mito'] < percent_mito]

    data.obs.loc[np.logical_and.reduce(filters), 'passed_qc'] = True

    var = data.var
    data = data[data.obs['passed_qc']]  # compute gene stats in space of filtered cells only

    var['n_cells'] = data.X.getnnz(axis=0)
    var['percent_cells'] = var['n_cells'] / data.shape[0]
    var['robust'] = var['percent_cells'] >= percent_cells
    var['highly_variable_features'] = var['robust']  # default all robust genes are "highly" variable


def get_filter_stats(data: 'AnnData') -> Tuple['pandas.DataFrame', 'pandas.DataFrame']:
    """
    TODO: documentation


    """

    # cell stats
    gb1 = data.obs.groupby('Channel')
    df_before = gb1.median()
    df_before = df_before.assign(total=gb1.size())
    df_before.rename(
        columns={
            'n_genes': 'median_n_genes_before',
            'n_counts': 'median_n_umis_before',
            'percent_mito': 'median_percent_mito_before',
        },
        inplace=True,
    )

    data = data[data.obs['passed_qc']] # focusing only on filtered cells

    gb2 = data.obs.groupby('Channel')
    df_after = gb2.median()
    df_after = df_after.assign(kept=gb2.size())
    df_after.rename(
        columns={
            'n_genes': 'median_n_genes',
            'n_counts': 'median_n_umis',
            'percent_mito': 'median_percent_mito',
        },
        inplace=True,
    )
    df_cells = pd.concat((df_before, df_after), axis=1, sort=False)
    df_cells.fillna(0, inplace=True)
    df_cells['kept'] = df_cells['kept'].astype(int)
    df_cells['filt'] = df_cells['total'] - df_cells['kept']
    df_cells = df_cells[
        [
            'kept',
            'median_n_genes',
            'median_n_umis',
            'median_percent_mito',
            'filt',
            'total',
            'median_n_genes_before',
            'median_n_umis_before',
            'median_percent_mito_before',
        ]
    ]
    df_cells.sort_values('kept', inplace=True)

    # gene stats
    idx = data.var['robust'] == False
    df_genes = pd.DataFrame(
        {
            'n_cells': data.var.loc[idx, 'n_cells'],
            'percent_cells': data.var.loc[idx, 'percent_cells'],
        }
    )
    df_genes.index.name = 'gene'
    df_genes.sort_values('n_cells', ascending=False, inplace=True)

    return df_cells, df_genes


def filter_data(data: 'AnnData') -> None:
    """ Filter data based on qc_metrics calculated
    TODO: Documentation
    """

    assert 'passed_qc' in data.obs
    data._inplace_subset_obs(data.obs['passed_qc'].values)
    data._inplace_subset_var((data.var['n_cells'] > 0).values)
    print(
        "After filteration, {nc} cells and {ng} genes are kept. Among {ng} genes, {nrb} genes are robust.".format(
            nc=data.shape[0], ng=data.shape[1], nrb=data.var["robust"].sum()
        )
    )


def generate_filter_plots(data: 'AnnData', plot_filt: str, plot_filt_figsize: str = None) -> None:
    """ This function generates filtration plots, only used in command line.
    """

    df_plot_before = data.obs[
        ['Channel', 'n_genes', 'n_counts', 'percent_mito']
    ].copy()
    df_plot_before.reset_index(drop=True, inplace=True)
    df_plot_before['status'] = 'original'

    data = data[data.obs['passed_qc']] # focusing only on filtered cells

    df_plot_after = data.obs[
        ['Channel', 'n_genes', 'n_counts', 'percent_mito']
    ].copy()
    df_plot_after.reset_index(drop=True, inplace=True)
    df_plot_after['status'] = 'filtered'
    df_plot = pd.concat((df_plot_before, df_plot_after), axis=0)
    
    from scCloud.plotting import plot_qc_violin
    figsize = None
    if plot_filt_figsize is not None:
        width, height = plot_filt_figsize.split(',')
        figsize = (int(width), int(height))

    plot_qc_violin(
        df_plot,
        'count',
        plot_filt + '.filt.UMI.pdf',
        xattr='Channel',
        hue='status',
        xlabel='Channel',
        split=True,
        linewidth=0,
        figsize=figsize,
    )

    plot_qc_violin(
        df_plot,
        'gene',
        plot_filt + '.filt.gene.pdf',
        xattr='Channel',
        hue='status',
        xlabel='Channel',
        split=True,
        linewidth=0,
        figsize=figsize,
    )

    plot_qc_violin(
        df_plot,
        'mito',
        plot_filt + '.filt.mito.pdf',
        xattr='Channel',
        hue='status',
        xlabel='Channel',
        split=True,
        linewidth=0,
        figsize=figsize,
    )

    print("Filtration plots are generated.")


def run_filter_data(
    data: 'AnnData',
    output_filt: str = None,
    plot_filt: str = None,
    plot_filt_figsize: Tuple[int, int] = None,
    mito_prefix: str = 'MT-',
    min_genes: int = 500,
    max_genes: int = 6000,
    min_umis: int = 100,
    max_umis: int = 600000,
    percent_mito: float = 0.1,
    percent_cells: float = 0.0005,
) -> None:
    """ This function is for command line use.
    TODO: Documentation
    """

    start = time.time()

    qc_metrics(data, mito_prefix, min_genes, max_genes, min_umis, max_umis, percent_mito, percent_cells)

    if output_filt is not None:
        writer = pd.ExcelWriter(output_filt + '.filt.xlsx', engine='xlsxwriter')
        df_cells, df_genes = get_filter_stats(data)
        df_cells.to_excel(writer, sheet_name='Cell filtration stats')
        df_genes.to_excel(writer, sheet_name='Gene filtration stats')
        writer.save()
        print("Filtration results are written.")

    if plot_filt is not None:
        generate_filter_plots(data, plot_filt, plot_filt_figsize)        

    filter_data(data)    

    end = time.time()
    print("filter_data is finished. Time spent = {:.2f}s.".format(end - start))


def log_norm(data: 'AnnData', norm_count: float = 1e5) -> None:
    """Normalization and then take log 
    TODO: Documentation
    """

    start = time.time()

    assert issparse(data.X)
    mat = data.X[:, data.var['robust'].values]
    scale = norm_count / mat.sum(axis=1).A1
    data.X.data *= np.repeat(scale, np.diff(data.X.indptr))
    data.X = data.X.log1p()

    end = time.time()
    print("Normalization is finished. Time spent = {:.2f}s.".format(end - start))


def select_features(data: 'AnnData', features: 'str' = None) -> str:
    """ Select a subset of features to form a new AnnData object with dense matrix. Store the matrix in data.uns with 'anndata_' prefix. Any key with 'anndata_' should be removed before writing out the disk.
    :param features: a keyword in data.var, which refers to a boolean array. None refers to all features

    TODO: Documentation
    """
    keyword = 'anndata_' + str(features)

    if keyword not in data.uns:    
        if features is not None:
            assert features in data.var
            data_dense = data[:, data.var[features]].copy()
        else:
            data_dense = data.copy()
        
        if issparse(data_dense.X):
            data_dense.X = data_dense.X.toarray()

        data.uns[keyword] = data_dense

    return keyword


def pca(data: 'AnnData', standardize: bool = True, max_value: float = 10, nPC: int = 50, random_state: int = 0, features: str = 'highly_variable_features') -> None:
    """Calculate PCA
    TODO: documentation. Feature selection time is not included.
    """

    keyword = select_features(data, features)
    
    start = time.time()

    data_dense = data.uns[keyword]

    if standardize:
        # scaler = StandardScaler(copy=False)
        # scaler.fit_transform(data_dense.X)
        X = data_dense.X 
        m1 = X.mean(axis = 0)
        psum = np.multiply(X, X).sum(axis = 0)
        std = ((psum - X.shape[0] * (m1 ** 2)) / (X.shape[0] - 1.0)) ** 0.5
        X -= m1
        X /= std

    if max_value is not None:
        data_dense.X[data_dense.X > max_value] = max_value
        data_dense.X[data_dense.X < -max_value] = -max_value

    pca = PCA(n_components=nPC, random_state=random_state)
    X_pca = pca.fit_transform(data_dense.X)

    data.obsm['X_pca'] = X_pca
    data.uns['PCs'] = pca.components_.T # cannot be varm because numbers of features are not the same
    data.uns['pca'] = {}
    data.uns['pca']['variance'] = pca.explained_variance_
    data.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_
    
    end = time.time()
    print("PCA is done. Time spent = {:.2f}s.".format(end - start))
