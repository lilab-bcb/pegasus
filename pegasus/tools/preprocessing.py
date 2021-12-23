import numpy as np
import pandas as pd

from scipy.sparse import issparse, csr_matrix
from threadpoolctl import threadpool_limits
from sklearn.decomposition import PCA, TruncatedSVD

from typing import List, Tuple, Union
from pegasusio import UnimodalData, MultimodalData, calc_qc_filters, DictWithDefault
from pegasus.tools import eff_n_jobs, normalize_by_count

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer, run_gc



def qc_metrics(
    data: MultimodalData,
    select_singlets: bool = False,
    remap_string: str = None,
    subset_string: str = None,
    min_genes: int = None,
    max_genes: int = None,
    min_umis: int = None,
    max_umis: int = None,
    mito_prefix: str = None,
    percent_mito: float = None,
) -> None:
    """Generate Quality Control (QC) metrics regarding cell barcodes on the dataset.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
       Use current selected modality in data, which should contain one RNA expression matrix.
    select_singlets: ``bool``, optional, default ``False``
        If select only singlets.
    remap_string: ``str``, optional, default ``None``
        Remap singlet names using <remap_string>, where <remap_string> takes the format "new_name_i:old_name_1,old_name_2;new_name_ii:old_name_3;...". For example, if we hashed 5 libraries from 3 samples sample1_lib1, sample1_lib2, sample2_lib1, sample2_lib2 and sample3, we can remap them to 3 samples using this string: "sample1:sample1_lib1,sample1_lib2;sample2:sample2_lib1,sample2_lib2". In this way, the new singlet names will be in metadata field with key 'assignment', while the old names will be kept in metadata field with key 'assignment.orig'.
    subset_string: ``str``, optional, default ``None``
        If select singlets, only select singlets in the <subset_string>, which takes the format "name1,name2,...". Note that if --remap-singlets is specified, subsetting happens after remapping. For example, we can only select singlets from sampe 1 and 3 using "sample1,sample3".
    min_genes: ``int``, optional, default: ``None``
       Only keep cells with at least ``min_genes`` genes.
    max_genes: ``int``, optional, default: ``None``
       Only keep cells with less than ``max_genes`` genes.
    min_umis: ``int``, optional, default: ``None``
       Only keep cells with at least ``min_umis`` UMIs.
    max_umis: ``int``, optional, default: ``None``
       Only keep cells with less than ``max_umis`` UMIs.
    mito_prefix: ``str``, optional, default: ``None``
       Prefix for mitochondrial genes.
    percent_mito: ``float``, optional, default: ``None``
       Only keep cells with percent mitochondrial genes less than ``percent_mito`` % of total counts.

    Returns
    -------
    ``None``

    Update ``data.obs``:

        * ``n_genes``: Total number of genes for each cell.
        * ``n_counts``: Total number of counts for each cell.
        * ``percent_mito``: Percent of mitochondrial genes for each cell.
        * ``passed_qc``: Boolean type indicating if a cell passes the QC process based on the QC metrics.
        * ``demux_type``: this column might be deleted if select_singlets is on.

    Examples
    --------
    >>> pg.qc_metrics(data, min_genes=500, max_genes=6000, mito_prefix="MT-", percent_mito=10)
    """
    if isinstance(data, MultimodalData):
        data = data.current_data()
    calc_qc_filters(data, select_singlets = select_singlets, remap_string = remap_string, subset_string = subset_string, min_genes = min_genes, max_genes = max_genes, min_umis = min_umis, max_umis = max_umis, mito_prefix = mito_prefix, percent_mito = percent_mito)


def get_filter_stats(data: MultimodalData, min_genes_before_filt: int = 100) -> pd.DataFrame:
    """Calculate filtration stats on cell barcodes.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Use current selected modality in data, which should contain one RNA expression matrix.
    min_genes_before_filt: ``int``, optional, default ``100``
        If raw data matrix is input, empty barcodes will dominate pre-filtration statistics. To avoid this, for raw matrix, only consider barcodes with at least <number> genes for pre-filtration condition.

    Returns
    -------
    df_cells: ``pandas.DataFrame``
        Data frame of stats on cell filtration.

    Examples
    --------
    >>> df = pg.get_filter_stats(data)
    """

    # cell stats
    if isinstance(data, MultimodalData):
        data = data.current_data()

    if "Channel" not in data.obs:
        data.obs["Channel"] = pd.Categorical([""] * data.shape[0])

    df = data.obs[data.obs["n_genes"] >= min_genes_before_filt] if data.obs["n_genes"].min() == 0 else data.obs
    gb1 = df.groupby("Channel")
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

    target_cols = np.array(["kept", "median_n_genes", "median_n_umis", "median_percent_mito", "filt", "total", "median_n_genes_before", "median_n_umis_before", "median_percent_mito_before"])
    target_cols = target_cols[np.isin(target_cols, df_cells.columns)]
    df_cells = df_cells[target_cols]
    df_cells.sort_values("kept", inplace=True)

    return df_cells


@timer(logger=logger)
def filter_data(data: MultimodalData, focus_list: List[str] = None) -> None:
    """ Filter data based on qc_metrics calculated in ``pg.qc_metrics``.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Use current selected modality in data, which should contain one RNA expression matrix.
    focus_list: ``List[str]``, optional, default None
        UnimodalData objects with keys in focus_list were qc_metrics marked. Filter them and make sure other modalities' barcodes are consistent with filtered barcodes. If focus_list is None and self._selected's modality is "rna", focus_list = [self._selected]

    Returns
    -------
    ``None``

    Update ``data`` with cells after filtration.

    Examples
    --------
    >>> pg.filter_data(data)
    """
    data.filter_data(focus_list = focus_list, cache_passqc = True, uns_white_list="norm_count,df_qcplot")


@timer(logger=logger)
def identify_robust_genes(data: MultimodalData, percent_cells: float = 0.05) -> None:
    """ Identify robust genes as candidates for HVG selection and remove genes that are not expressed in any cells.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Use current selected modality in data, which should contain one RNA expression matrix.
    percent_cells: ``float``, optional, default: ``0.05``
       Only assign genes to be ``robust`` that are expressed in at least ``percent_cells`` % of cells.

    Returns
    -------
    ``None``

    Update ``data.var``:

        * ``n_cells``: Total number of cells in which each gene is measured.
        * ``percent_cells``: Percent of cells in which each gene is measured.
        * ``robust``: Boolean type indicating if a gene is robust based on the QC metrics.
        * ``highly_variable_features``: Boolean type indicating if a gene is a highly variable feature. By default, set all robust genes as highly variable features.

    Examples
    --------
    >>> pg.identify_robust_genes(data, percent_cells = 0.05)
    """
    if isinstance(data, MultimodalData):
        data = data.current_data()

    prior_n = data.shape[1]

    if issparse(data.X):
        data.var["n_cells"] = data.X.getnnz(axis=0)
        data._inplace_subset_var(data.var["n_cells"] > 0)
        data.var["percent_cells"] = (data.var["n_cells"] / data.shape[0]) * 100
        data.var["robust"] = data.var["percent_cells"] >= percent_cells
    else:
        data.var["robust"] = True

    data.var["highly_variable_features"] = data.var["robust"]  # default all robust genes are "highly" variable
    logger.info(f"After filtration, {data.shape[1]}/{prior_n} genes are kept. Among {data.shape[1]} genes, {data.var['robust'].sum()} genes are robust.")


def _generate_filter_plots(
    unidata: UnimodalData, plot_filt: str, plot_filt_figsize: str = None, min_genes_before_filt: int = 100
) -> None:
    """ This function generates filtration plots, only used in command line.
    """
    group_key = unidata.get_uid()

    from pegasus.plotting import qcviolin

    kwargs = {"return_fig": True, "dpi": 500}
    if plot_filt_figsize is not None:
        width, height = plot_filt_figsize.split(",")
        kwargs["panel_size"] = (int(width), int(height))

    fig = qcviolin(unidata, "count", **kwargs)
    fig.savefig(f"{plot_filt}.{group_key}.filt.UMI.pdf")

    fig = qcviolin(unidata, "gene", **kwargs)
    fig.savefig(f"{plot_filt}.{group_key}.filt.gene.pdf")

    fig = qcviolin(unidata, "mito", **kwargs)
    if fig is not None:
        fig.savefig(f"{plot_filt}.{group_key}.filt.mito.pdf")

    logger.info("Filtration plots are generated.")


@timer(logger=logger)
def _run_filter_data(
    data: MultimodalData,
    focus_list: List[str] = None,
    output_filt: str = None,
    plot_filt: str = None,
    plot_filt_figsize: Tuple[int, int] = None,
    min_genes_before_filt: int = 100,
    select_singlets: bool = False,
    remap_string: str = None,
    subset_string: str = None,
    min_genes: int = 500,
    max_genes: int = 6000,
    min_umis: int = None,
    max_umis: int = None,
    mito_prefix: str = "MT-",
    percent_mito: float = 20.0,
    percent_cells: float = 0.05,
) -> None:
    """ This function is for command line use.
    """
    if focus_list is None:
        focus_list = [data.current_key()]

    mito_dict = DictWithDefault(mito_prefix)
    for key in focus_list:
        unidata = data.get_data(key)

        qc_metrics(
            unidata,
            select_singlets,
            remap_string,
            subset_string,
            min_genes,
            max_genes,
            min_umis,
            max_umis,
            mito_dict.get(unidata.get_genome()),
            percent_mito,
        )

        if output_filt is not None:
            group_key = unidata.get_uid()
            writer = pd.ExcelWriter(f"{output_filt}.{group_key}.filt.xlsx", engine="xlsxwriter")
            df_cells = get_filter_stats(unidata, min_genes_before_filt = min_genes_before_filt)
            df_cells.to_excel(writer, sheet_name="Cell filtration stats")
            writer.save()
            logger.info(f"Filtration results for {group_key} are written.")

        if plot_filt is not None:
            _generate_filter_plots(unidata, plot_filt, plot_filt_figsize = plot_filt_figsize, min_genes_before_filt = min_genes_before_filt)

    filter_data(data, focus_list = focus_list)

    for key in focus_list:
        unidata = data.get_data(key)
        identify_robust_genes(unidata, percent_cells = percent_cells)


@timer(logger=logger)
def log_norm(
    data: Union[MultimodalData, UnimodalData],
    norm_count: float = 1e5, 
    backup_matrix: str = "raw.X",
) -> None:
    """Normalization, and then apply natural logarithm to the data.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Use current selected modality in data, which should contain one RNA expression matrix.

    norm_count: ``int``, optional, default: ``1e5``.
        Total counts of one cell after normalization.

    backup_matrix: ``str``, optional, default: ``raw.X``.
        The key name of the backup count matrix, usually the raw counts.

    Returns
    -------
    ``None``

    Update ``data.X`` with count matrix after log-normalization. In addition, back up the original count matrix as ``backup_matrix``.

    In case of rerunning normalization while ``backup_matrix`` already exists, use ``backup_matrix`` instead of ``data.X`` for normalization.

    Examples
    --------
    >>> pg.log_norm(data)
    """
    if isinstance(data, MultimodalData):
        data = data.current_data()

    assert data.get_modality() in  {"rna", "visium"}

    if backup_matrix not in data.list_keys():
        data.add_matrix(backup_matrix, data.X)
        data.X = data.X.astype(np.float32)  # force copy
    else:
        # The case of rerunning log_norm. Use backup matrix as source.
        data.X = data.get_matrix(backup_matrix).astype(np.float32)  # force copy
        logger.warning("Rerun log-normalization. Use the raw counts in backup instead.")

    data.obs["scale"] = normalize_by_count(data.X, data.var["robust"].values, norm_count, True)
    data.uns["norm_count"] = norm_count


@run_gc
def select_features(
    data: Union[MultimodalData, UnimodalData],
    features: str = "highly_variable_features",
    standardize: bool = True,
    max_value: float = 10.0,
) -> str:
    """ Subset the features and store the resulting matrix in dense format in data.uns with `'_tmp_fmat_'` prefix, with the option of standardization and truncating based on max_value. `'_tmp_fmat_*'` will be removed before writing out the disk.
    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    features: ``str``, optional, default: ``highly_variable_features``.
        a keyword in ``data.var``, which refers to a boolean array. If ``None``, all features will be selected.
    standardize: ``bool``, optional, default: ``True``.
        Whether to scale the data to unit variance and zero mean.
    max_value: ``float``, optional, default: ``10``.
        The threshold to truncate data after scaling. If ``None``, do not truncate.
    Returns
    -------
    keyword: ``str``
        The keyword in ``data.uns`` referring to the features selected.
    Update ``data.uns`` if needed:
        * ``data.uns[keyword]``: A submatrix of the data containing features selected.
    Examples
    --------
    >>> pg.select_features(data)
    """
    keyword = "_tmp_fmat_" + str(features)  # fmat: feature matrix
    batch_ls_key = "_tmp_ls_" + str(features)

    if batch_ls_key in data.uns:
        X = data.uns[keyword]
    else:
        if features is not None:
            assert features in data.var
            X = data.X[:, data.var[features].values]
        else:
            X = data.X
        
        from pegasus.tools import slicing
        X = slicing(X, copy=True)
        
    if standardize:
        m1 = X.mean(axis=0)
        psum = np.multiply(X, X).sum(axis=0)
        std = ((psum - X.shape[0] * (m1 ** 2)) / (X.shape[0] - 1.0)) ** 0.5
        std[std == 0] = 1
        X -= m1
        X /= std
        data.uns["stdzn_mean"] = m1
        data.uns["stdzn_std"] = std
        data.uns.pop(batch_ls_key, None)

    if max_value is not None:
        data.uns["stdzn_max_value"] = max_value
        X[X > max_value] = max_value
        X[X < -max_value] = -max_value
        data.uns.pop(batch_ls_key, None)

    data.uns[keyword] = X

    return keyword


@timer(logger=logger)
def pca(
    data: Union[MultimodalData, UnimodalData],
    n_components: int = 50,
    features: str = "highly_variable_features",
    standardize: bool = True,
    max_value: float = 10.0,
    n_jobs: int = -1,
    random_state: int = 0,
) -> None:
    """Perform Principle Component Analysis (PCA) to the data.

    The calculation uses *scikit-learn* implementation.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``50``.
        Number of Principal Components to get.

    features: ``str``, optional, default: ``"highly_variable_features"``.
        Keyword in ``data.var`` to specify features used for PCA.

    standardize: ``bool``, optional, default: ``True``.
        Whether to scale the data to unit variance and zero mean.

    max_value: ``float``, optional, default: ``10``.
        The threshold to truncate data after scaling. If ``None``, do not truncate.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

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

        * ``data.uns["pca_features"]``: Record the features used to generate PCA components.

    Examples
    --------
    >>> pg.pca(data)
    """
    keyword = select_features(data, features=features, standardize=standardize, max_value=max_value)
    X = data.uns[keyword].astype(np.float64) # float64 to avoid precision issues and make results more reproducible across platforms
    pca = PCA(n_components=n_components, random_state=random_state) # use auto solver, default is randomized for large datasets

    n_jobs = eff_n_jobs(n_jobs)
    with threadpool_limits(limits = n_jobs):
        X_pca = pca.fit_transform(X)

    data.obsm["X_pca"] = np.ascontiguousarray(X_pca, dtype=np.float32)
    data.uns["PCs"] = np.ascontiguousarray(pca.components_.T, dtype=np.float32)  # cannot be varm because numbers of features are not the same
    data.uns["pca"] = {}
    data.uns["pca"]["variance"] = pca.explained_variance_
    data.uns["pca"]["variance_ratio"] = pca.explained_variance_ratio_
    data.uns["pca_features"] = features # record which feature to use


@timer(logger=logger)
def pc_transform(
    data: Union[MultimodalData, UnimodalData],
    X: Union[csr_matrix, np.ndarray]
) -> np.ndarray:
    """Transform a new count matrix X using existing PCs contained in data.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with PCs calculated

    X: ``Union[csr_matrix, np.ndarray]``
        New input count matrix to be transformed. Number of genes must be the same as in data

    Returns
    -------
    Transformed PC coordinates for X.

    Examples
    --------
    >>> pg.pc_transform(data, X_new)
    """
    assert data.shape[1] == X.shape[1]
    if ("norm_count" not in data.uns) or ("robust" not in data.var) or ("pca_features" not in data.uns) or ("PCs" not in data.uns):
        raise ValueError("pc_transform cannot be applied before identify_robust_gene, highly_variable_features and pca are called!")

    # normalize and log transform X
    X = X.astype(np.float32) # convert to float
    _ = normalize_by_count(X, data.var["robust"].values, data.uns["norm_count"], True)

    # select features
    if data.uns["pca_features"] is not None:
        X = X[:, data.var[data.uns["pca_features"]].values]
    from pegasus.tools import slicing
    X = slicing(X, copy=True)

    # standardization
    if "stdzn_mean" in data.uns:
        assert "stdzn_std" in data.uns
        X -= data.uns["stdzn_mean"]
        X /= data.uns["stdzn_std"]
    if "stdzn_max_value" in data.uns:
        max_value = data.uns["stdzn_max_value"]
        X[X > max_value] = max_value
        X[X < -max_value] = -max_value

    # transform to PC coordinates
    X_pca = np.matmul(X, data.uns["PCs"])
    return X_pca


@timer(logger=logger)
def regress_out(
    data: Union[MultimodalData, UnimodalData],
    attrs: List[str],
    rep: str = 'pca',
) -> str:
    """Regress out effects due to specific observational attributes.

    Parameters
    ----------
    data: ``MultimodalData`` or ``UnimodalData`` object
        Annotated data matrix with rows for cells and columns for genes.

    attrs: ``List[str]``
        List of numeric cell attributes to be regressed out. They must exist in ``data.obs`` field.

    rep: ``str``, optional, default: ``pca``
        This is to specify which embedding to be used for regressing out.
        The key ``'X_'+rep`` must exist in ``data.obsm`` field. By default, use PCA embedding.

    Returns
    -------
    res_key: ``str``
        The key to the resulting new embedding matrix in ``data.obsm``. It's ``'X_'+rep+'_regressed'``.

    Update ``data.obsm``:
        * ``data.obsm[pca_key]``: The PCA matrix with effects of attributes specified regressed out.

    Examples
    --------
    >>> pg.regress_out(data, attrs=['G1/S', 'G2/M'])
    """
    n_components = data.obsm[f'X_{rep}'].shape[1]

    from pandas.api.types import is_numeric_dtype
    for attr in attrs:
        if not is_numeric_dtype(data.obs[attr]):
            raise TypeError(f"Cell attribute '{attr}' is not numeric. For regressing out, All attributes used must be numeric!")

    X = data.obs[attrs]

    from sklearn.linear_model import LinearRegression

    response_list = []
    for i in range(n_components):
        pc = data.obsm[f'X_{rep}'][:, i]
        model = LinearRegression().fit(X, pc)
        y_pred = model.predict(X)
        resid = pc - y_pred
        response_list.append(resid)

    pca_key = f'{rep}_regressed'
    data.obsm[f'X_{pca_key}'] = np.ascontiguousarray(np.vstack(response_list).T.astype(data.obsm[f'X_{rep}'].dtype), dtype=np.float32)

    return pca_key


@timer(logger=logger)
def tsvd(
    data: MultimodalData,
    n_components: int = 51,
    features: str = "robust",
    n_jobs: int = -1,
    random_state: int = 0,
) -> None:
    """Perform Truncated Singular Value Decomposition (TSVD) to the data.

    The calculation uses *scikit-learn* implementation.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.

    n_components: ``int``, optional, default: ``51``.
        Number of components to compute. Because we do not remove means, the first component might be related to the mean values and thus we ask for one more component than PCA as default (i.e. 51).

    features: ``str``, optional, default: ``"robust"``.
        Keyword in ``data.var`` to specify features used for TSVD.

    robust: ``bool``, optional, default: ``False``.
        If true, use 'arpack' instead of 'randomized'.

    n_jobs : `int`, optional (default: -1)
        Number of threads to use. -1 refers to using all physical CPU cores.

    random_state: ``int``, optional, default: ``0``.
        Random seed to be set for reproducing result.

    Returns
    -------
    ``None``.

    Update ``data.obsm``:

        * ``data.obsm["X_tsvd"]``: TSVD coordinates of the data.

    Update ``data.uns``:

        * ``data.uns["TSVDs"]``: The TSVD components containing the loadings.

    Examples
    --------
    >>> pg.tsvd(data)
    """
    X = data.X
    if (features is not None) and (data.var[features].sum() < data.shape[1]):
        X = X[:, data.var[features].values]
    X = X.astype(np.float64) # for reproducible purpose

    tsvd = TruncatedSVD(n_components = n_components, random_state = random_state)

    n_jobs = eff_n_jobs(n_jobs)
    with threadpool_limits(limits = n_jobs):
        X_tsvd = tsvd.fit_transform(X)

    data.obsm["X_tsvd"] = np.ascontiguousarray(X_tsvd, dtype=np.float32)
    data.uns["TSVDs"] = np.ascontiguousarray(tsvd.components_.T, dtype=np.float32)  # cannot be varm because numbers of features are not the same
    data.uns["tsvd"] = {}
    data.uns["tsvd"]["variance"] = tsvd.explained_variance_
    data.uns["tsvd"]["variance_ratio"] = tsvd.explained_variance_ratio_
    data.uns["tsvd_features"] = features # record which feature to use


@timer(logger=logger)
def tsvd_transform(
    data: MultimodalData,
    X: Union[csr_matrix, np.ndarray]
) -> np.ndarray:
    """Transform a new count matrix X using existing TSVD components contained in data.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with PCs calculated

    X: ``Union[csr_matrix, np.ndarray]``
        New input count matrix to be transformed. Number of genes must be the same as in data

    Returns
    -------
    Transformed TSVD coordinates for X.

    Examples
    --------
    >>> pg.tsvd_transform(data, X_new)
    """
    assert data.shape[1] == X.shape[1]
    if ("norm_count" not in data.uns) or ("robust" not in data.var) or ("tsvd_features" not in data.uns) or ("TSVDs" not in data.uns):
        raise ValueError("tsvd_transform cannot be applied before identify_robust_gene, and tsvd are called!")

    # normalize and log transform X
    X = X.astype(np.float32) # convert to float
    _ = normalize_by_count(X, data.var["robust"].values, data.uns["norm_count"], True)

    # select features
    if (data.uns["tsvd_features"] is not None) and (data.var[data.uns["tsvd_features"]].sum() < X.shape[1]):
        X = X[:, data.var[data.uns["tsvd_features"]].values]

    # transform to TSVD coordinates
    from sklearn.utils.extmath import safe_sparse_dot
    X_tsvd = safe_sparse_dot(X, data.uns["TSVDs"])
    return X_tsvd
