import time
import numpy as np
import pandas as pd

from scipy.sparse import issparse

from sklearn.decomposition import PCA

from typing import Tuple
from anndata import AnnData

import logging

logger = logging.getLogger("pegasus")
from pegasus.utils import decorators as pg_deco


@pg_deco.GCCollect()
def qc_metrics(
    data: AnnData,
    mito_prefix: str = "MT-",
    min_genes: int = 500,
    max_genes: int = 6000,
    min_umis: int = 100,
    max_umis: int = 600000,
    percent_mito: float = 10.0,
    percent_cells: float = 0.05,
) -> None:
    """Generate Quality Control (QC) metrics on the dataset.

    Parameters
    ----------
    data: ``anndata.AnnData``
       Annotated data matrix with rows for cells and columns for genes.
    mito_prefix: ``str``, optional, default: ``"MT-"``
       Prefix for mitochondrial genes.
    min_genes: ``int``, optional, default: ``500``
       Only keep cells with at least ``min_genes`` genes.
    max_genes: ``int``, optional, default: ``6000``
       Only keep cells with less than ``max_genes`` genes.
    min_umis: ``int``, optional, default: ``100``
       Only keep cells with at least ``min_umis`` UMIs.
    max_umis: ``int``, optional, default: ``600,000``
       Only keep cells with less than ``max_umis`` UMIs.
    percent_mito: ``float``, optional, default: ``10.0``
       Only keep cells with percent mitochondrial genes less than ``percent_mito`` % of total counts.
    percent_cells: ``float``, optional, default: ``0.05``
       Only assign genes to be ``robust`` that are expressed in at least ``percent_cells`` % of cells.

    Returns
    -------
    ``None``

    Update ``data.obs``:

        * ``n_genes``: Total number of genes for each cell.
        * ``n_counts``: Total number of counts for each cell.
        * ``percent_mito``: Percent of mitochondrial genes for each cell.
        * ``passed_qc``: Boolean type indicating if a cell passes the QC process based on the QC metrics.

    Update ``data.var``:

        * ``n_cells``: Total number of cells in which each gene is measured.
        * ``percent_cells``: Percent of cells in which each gene is measured.
        * ``robust``: Boolean type indicating if a gene is robust based on the QC metrics.
        * ``highly_variable_features``: Boolean type indicating if a gene is a highly variable feature. By default, set all robust genes as highly variable features.

    Examples
    --------
    >>> pg.qcmetrics(adata)
    """

    data.obs["passed_qc"] = False

    data.obs["n_genes"] = data.X.getnnz(axis=1)
    data.obs["n_counts"] = data.X.sum(axis=1).A1

    mito_prefixes = mito_prefix.split(",")

    def startswith(name):
        for prefix in mito_prefixes:
            if name.startswith(prefix):
                return True
        return False

    mito_genes = data.var_names.map(startswith).values.nonzero()[0]
    data.obs["percent_mito"] = (
        data.X[:, mito_genes].sum(axis=1).A1
        / np.maximum(data.obs["n_counts"].values, 1.0)
    ) * 100

    # Assign passed_qc
    filters = [
        data.obs["n_genes"] >= min_genes,
        data.obs["n_genes"] < max_genes,
        data.obs["n_counts"] >= min_umis,
        data.obs["n_counts"] < max_umis,
        data.obs["percent_mito"] < percent_mito,
    ]

    data.obs.loc[np.logical_and.reduce(filters), "passed_qc"] = True

    var = data.var
    data = data[
        data.obs["passed_qc"]
    ]  # compute gene stats in space of filtered cells only

    var["n_cells"] = data.X.getnnz(axis=0)
    var["percent_cells"] = (var["n_cells"] / data.shape[0]) * 100
    var["robust"] = var["percent_cells"] >= percent_cells
    var["highly_variable_features"] = var[
        "robust"
    ]  # default all robust genes are "highly" variable


def get_filter_stats(data: AnnData) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate filtration stats on cell barcodes and genes, respectively.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    Returns
    -------
    df_cells: ``pandas.DataFrame``
        Data frame of stats on cell filtration.

    df_genes: ``pandas.DataFrame``
        Data frame of stats on gene filtration.

    Examples
    --------
    >>> pg.get_filter_stats(adata)
    """

    # cell stats
    gb1 = data.obs.groupby("Channel")
    df_before = gb1.median()
    df_before = df_before.assign(total=gb1.size())
    df_before.rename(
        columns={
            "n_genes": "median_n_genes_before",
            "n_counts": "median_n_umis_before",
            "percent_mito": "median_percent_mito_before",
        },
        inplace=True,
    )

    # focusing only on filtered cells
    gb2 = data.obs[data.obs["passed_qc"]].groupby("Channel")
    df_after = gb2.median()
    df_after = df_after.assign(kept=gb2.size())
    df_after.rename(
        columns={
            "n_genes": "median_n_genes",
            "n_counts": "median_n_umis",
            "percent_mito": "median_percent_mito",
        },
        inplace=True,
    )
    df_cells = pd.concat((df_before, df_after), axis=1, sort=False)
    df_cells.fillna(0, inplace=True)
    df_cells["kept"] = df_cells["kept"].astype(int)
    df_cells["filt"] = df_cells["total"] - df_cells["kept"]
    df_cells = df_cells[
        [
            "kept",
            "median_n_genes",
            "median_n_umis",
            "median_percent_mito",
            "filt",
            "total",
            "median_n_genes_before",
            "median_n_umis_before",
            "median_percent_mito_before",
        ]
    ]
    df_cells.sort_values("kept", inplace=True)

    # gene stats
    idx = data.var["robust"] == False
    df_genes = pd.DataFrame(
        {
            "n_cells": data.var.loc[idx, "n_cells"],
            "percent_cells": data.var.loc[idx, "percent_cells"],
        }
    )
    df_genes.index.name = "gene"
    df_genes.sort_values("n_cells", ascending=False, inplace=True)

    return df_cells, df_genes


def filter_data(data: AnnData) -> None:
    """ Filter data based on qc_metrics calculated in ``pg.qc_metrics``.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    Returns
    -------
    ``None``

    Update ``data`` with cells and genes after filtration.

    Examples
    --------
    >>> pg.filter_data(adata)
    """

    assert "passed_qc" in data.obs
    prior_shape = data.shape
    data._inplace_subset_obs(data.obs["passed_qc"].values)
    data._inplace_subset_var((data.var["n_cells"] > 0).values)
    logger.info(
        "After filtration, {nc}/{ncp} cells and {ng}/{ngp} genes are kept. Among {ng} genes, {nrb} genes are robust.".format(
            nc=data.shape[0],
            ng=data.shape[1],
            ncp=prior_shape[0],
            ngp=prior_shape[1],
            nrb=data.var["robust"].sum(),
        )
    )


def generate_filter_plots(
    data: AnnData, plot_filt: str, plot_filt_figsize: str = None
) -> None:
    """ This function generates filtration plots, only used in command line.
    """

    df_plot_before = data.obs[["Channel", "n_genes", "n_counts", "percent_mito"]].copy()
    df_plot_before.reset_index(drop=True, inplace=True)
    df_plot_before["status"] = "original"

    data = data[data.obs["passed_qc"]]  # focusing only on filtered cells

    df_plot_after = data.obs[["Channel", "n_genes", "n_counts", "percent_mito"]].copy()
    df_plot_after.reset_index(drop=True, inplace=True)
    df_plot_after["status"] = "filtered"
    df_plot = pd.concat((df_plot_before, df_plot_after), axis=0)

    from pegasus.plotting import plot_qc_violin

    figsize = None
    if plot_filt_figsize is not None:
        width, height = plot_filt_figsize.split(",")
        figsize = (int(width), int(height))

    plot_qc_violin(
        df_plot,
        "count",
        plot_filt + ".filt.UMI.pdf",
        xattr="Channel",
        hue="status",
        xlabel="Channel",
        split=True,
        linewidth=0,
        figsize=figsize,
    )

    plot_qc_violin(
        df_plot,
        "gene",
        plot_filt + ".filt.gene.pdf",
        xattr="Channel",
        hue="status",
        xlabel="Channel",
        split=True,
        linewidth=0,
        figsize=figsize,
    )

    plot_qc_violin(
        df_plot,
        "mito",
        plot_filt + ".filt.mito.pdf",
        xattr="Channel",
        hue="status",
        xlabel="Channel",
        split=True,
        linewidth=0,
        figsize=figsize,
    )

    logger.info("Filtration plots are generated.")


@pg_deco.TimeLogger()
def run_filter_data(
    data: AnnData,
    output_filt: str = None,
    plot_filt: str = None,
    plot_filt_figsize: Tuple[int, int] = None,
    mito_prefix: str = "MT-",
    min_genes: int = 500,
    max_genes: int = 6000,
    min_umis: int = 100,
    max_umis: int = 600000,
    percent_mito: float = 10.0,
    percent_cells: float = 0.05,
) -> None:
    """ This function is for command line use.
    """

    qc_metrics(
        data,
        mito_prefix,
        min_genes,
        max_genes,
        min_umis,
        max_umis,
        percent_mito,
        percent_cells,
    )

    if output_filt is not None:
        writer = pd.ExcelWriter(output_filt + ".filt.xlsx", engine="xlsxwriter")
        df_cells, df_genes = get_filter_stats(data)
        df_cells.to_excel(writer, sheet_name="Cell filtration stats")
        df_genes.to_excel(writer, sheet_name="Gene filtration stats")
        writer.save()
        logger.info("Filtration results are written.")

    if plot_filt is not None:
        generate_filter_plots(data, plot_filt, plot_filt_figsize)

    filter_data(data)


@pg_deco.TimeLogger()
@pg_deco.GCCollect()
def log_norm(data: AnnData, norm_count: float = 1e5) -> None:
    """Normalization, and then apply natural logarithm to the data.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    norm_count: ``int``, optional, default: ``1e5``.
        Total count of cells after normalization.

    Returns
    -------
    ``None``

    Update ``data.X`` with count matrix after log-normalization.

    Examples
    --------
    >>> pg.log_norm(adata)
    """

    assert issparse(data.X)
    mat = data.X[:, data.var["robust"].values]
    scale = norm_count / mat.sum(axis=1).A1
    data.X.data *= np.repeat(scale, np.diff(data.X.indptr))
    data.X.data = np.log1p(data.X.data) # faster than data.X.log1p()


def select_features(data: AnnData, features: str = None) -> str:
    """ Subset the features and store the resulting matrix in dense format in data.uns with `'fmat_'` prefix. `'fmat_*'` will be removed before writing out the disk.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    features: ``str``, optional, default: ``None``.
        a keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.

    Returns
    -------
    keyword: ``str``
        The keyword in ``data.uns`` referring to the features selected.

    Update ``data.uns``:

        * ``data.uns[keyword]``: A submatrix of the data containing features selected.

    Examples
    --------
    >>> pg.select_features(adata)
    """
    keyword = "fmat_" + str(features)  # fmat: feature matrix

    if keyword not in data.uns:
        if features is not None:
            assert features in data.var
            fmat = data.X[:, data.var[features].values]
        else:
            fmat = data.X

        if issparse(fmat):
            data.uns[keyword] = fmat.toarray()
        else:
            data.uns[keyword] = fmat.copy()

    return keyword


@pg_deco.GCCollect()
def pca(
    data: AnnData,
    n_components: int = 50,
    features: str = "highly_variable_features",
    standardize: bool = True,
    max_value: float = 10,
    robust: bool = False,
    random_state: int = 0,
) -> None:
    """Perform Principle Component Analysis (PCA) to the data.

    The calculation uses *scikit-learn* implementation.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``50``.
        Number of Principal Components to get.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for PCA.

    standardize: ``bool``, optional, default: ``True``.
        Whether to scale the data to unit variance and zero mean or not.

    max_value: ``float``, optional, default: ``10``.
        The threshold to truncate data after scaling. If ``None``, do not truncate.

    robust: ``bool``, optional, default: ``False``.
        If true, use 'arpack' instead of 'randomized' for large sparse matrices (i.e. max(X.shape) > 500 and n_components < 0.8 * min(X.shape))

    random_state: ``int``, optional, default: ``0``.
        Random seed to be set for reproducing result.


    Returns
    -------
    ``None``.

    Update ``data.obsm``:

        * ``data.obsm["X_pca"]``: PCA matrix of the data.

    Update ``data.uns``:

        * ``data.uns["PCs"]``: The principal components containing the loadings.

        * ``data.uns["pca_variance"]``: Explained variance, i.e. the eigenvalues of the covariance matrix.

        * ``data.uns["pca_variance_ratio"]``: Ratio of explained variance.

    Examples
    --------
    >>> pg.pca(adata)
    """

    keyword = select_features(data, features)

    start = time.perf_counter()

    X = data.uns[keyword]

    if standardize:
        # scaler = StandardScaler(copy=False)
        # scaler.fit_transform(X)
        m1 = X.mean(axis=0)
        psum = np.multiply(X, X).sum(axis=0)
        std = ((psum - X.shape[0] * (m1 ** 2)) / (X.shape[0] - 1.0)) ** 0.5
        std[std == 0] = 1
        X -= m1
        X /= std

    if max_value is not None:
        X[X > max_value] = max_value
        X[X < -max_value] = -max_value

    pca = PCA(n_components=n_components, random_state=random_state)
    if robust:
        svd_solver = "arpack" if max(X.shape) > 500 and n_components < 0.8 * min(X.shape) else "full"
        pca = PCA(n_components=n_components, random_state=random_state, svd_solver=svd_solver)

    X_pca = pca.fit_transform(X)

    data.obsm["X_pca"] = X_pca
    data.uns[
        "PCs"
    ] = pca.components_.T  # cannot be varm because numbers of features are not the same
    data.uns["pca"] = {}
    data.uns["pca"]["variance"] = pca.explained_variance_
    data.uns["pca"]["variance_ratio"] = pca.explained_variance_ratio_

    end = time.perf_counter()
    logger.info("PCA is done. Time spent = {:.2f}s.".format(end - start))
