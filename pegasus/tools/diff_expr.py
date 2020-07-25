import time
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix, csc_matrix
from joblib import Parallel, delayed, effective_n_jobs
from statsmodels.stats.multitest import fdrcorrection as fdr
from collections import defaultdict

from typing import List, Tuple, Dict, Union

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer, MultimodalData, UnimodalData


def calc_basic_stat(
    clust_id: str,
    data: List[float],
    indices: List[int],
    indptr: List[int],
    shape: Tuple[int, int],
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    sum_vec: List[float],
    cnt_vec: List[int],
    verbose: bool,
) -> pd.DataFrame:
    """ Calcualte basic statistics for one cluster
    """
    zero_vec = np.zeros(shape[1], dtype=np.float32)

    # recover sparse matrix
    mat = csr_matrix((data, indices, indptr), shape=shape)
    mask = cluster_labels == clust_id
    mat_clust = mat[mask]

    if cond_labels is None:
        n1 = mat_clust.shape[0]
        n2 = shape[0] - n1

        sum_clu = mat_clust.sum(axis=0).A1
        mean1 = sum_clu / n1 if n1 > 0 else zero_vec
        mean2 = (sum_vec - sum_clu) / n2 if n2 > 0 else zero_vec
        mean2[mean2 < 0.0] = 0.0

        nonzeros = mat_clust.getnnz(axis=0)
        percents = (nonzeros / n1 * 100.0).astype(np.float32) if n1 > 0 else zero_vec
        percents_other = (
            ((cnt_vec - nonzeros) / n2 * 100.0).astype(np.float32)
            if n2 > 0
            else zero_vec
        )
    else:
        cond1 = cond_labels.categories[0]
        cond_labs = cond_labels[mask]
        mask2 = cond_labs == cond1

        mat_cond1 = mat_clust[mask2]
        mat_cond2 = mat_clust[~mask2]
        n1 = mat_cond1.shape[0]
        n2 = mat_cond2.shape[0]

        mean1 = mat_cond1.mean(axis=0).A1
        mean2 = mat_cond2.mean(axis=0).A1

        percents = (
            (mat_cond1.getnnz(axis=0) / n1 * 100.0).astype(np.float32)
            if n1 > 0
            else zero_vec
        )
        percents_other = (
            (mat_cond2.getnnz(axis=0) / n2 * 100.0).astype(np.float32)
            if n2 > 0
            else zero_vec
        )

    # calculate log_fold_change and WAD, Weighted Average Difference, https://almob.biomedcentral.com/articles/10.1186/1748-7188-3-8
    log_fold_change = mean1 - mean2
    x_avg = (mean1 + mean2) / 2
    x_max = x_avg.max()
    x_min = x_avg.min() - 0.001  # to avoid divide by zero
    weights = (x_avg - x_min) / (x_max - x_min)
    WADs = log_fold_change * weights
    # calculate percent fold change
    idx = percents > 0.0
    idx_other = percents_other > 0.0
    percent_fold_change = np.zeros(shape[1], dtype=np.float32)
    percent_fold_change[(~idx) & (~idx_other)] = np.nan
    percent_fold_change[idx & (~idx_other)] = np.inf
    percent_fold_change[idx_other] = percents[idx_other] / percents_other[idx_other]

    df = pd.DataFrame(
        {
            "mean_logExpr:{0}".format(clust_id): mean1,
            "mean_logExpr_other:{0}".format(clust_id): mean2,
            "log_fold_change:{0}".format(clust_id): log_fold_change,
            "fold_change:{0}".format(clust_id): np.exp(log_fold_change),
            "percentage:{0}".format(clust_id): percents,
            "percentage_other:{0}".format(clust_id): percents_other,
            "percentage_fold_change:{0}".format(clust_id): percent_fold_change,
            "WAD_score:{0}".format(clust_id): WADs,
        },
        index=gene_names,
    )

    if verbose:
        logger.info("calc_basic_stat finished for cluster {0}.".format(clust_id))

    return df


def collect_basic_statistics(
    X: csr_matrix,
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    n_jobs: int,
    temp_folder: str,
    verbose: bool,
) -> List[pd.DataFrame]:
    """ Collect basic statistics, triggering calc_basic_stat in parallel
    """
    start = time.perf_counter()

    sum_vec = cnt_vec = None
    if cond_labels is None:
        sum_vec = X.sum(axis=0).A1
        cnt_vec = X.getnnz(axis=0)

    result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
        delayed(calc_basic_stat)(
            clust_id,
            X.data,
            X.indices,
            X.indptr,
            X.shape,
            cluster_labels,
            cond_labels,
            gene_names,
            sum_vec,
            cnt_vec,
            verbose,
        )
        for clust_id in cluster_labels.categories
    )

    end = time.perf_counter()
    if verbose:
        logger.info(
            "Collecting basic statistics is done. Time spent = {:.2f}s.".format(
                end - start
            )
        )

    return result_list


def calc_auc(
    clust_id: str,
    data: List[float],
    indices: List[int],
    indptr: List[int],
    shape: Tuple[int, int],
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    verbose: bool,
) -> pd.DataFrame:
    """ Calculate AUROC for one cluster
    """
    import sklearn.metrics as sm

    csc_mat = csc_matrix((data, indices, indptr), shape=shape)
    mask = cluster_labels == clust_id

    auroc = np.zeros(shape[1], dtype=np.float32)
    # aupr = np.zeros(shape[1], dtype = np.float32)

    if cond_labels is None:
        exprs = np.zeros(shape[0])
        y_true = mask
        n1 = mask.sum()
        n2 = shape[0] - n1
    else:
        exprs = None
        cond1 = cond_labels.categories[0]
        cond_labs = cond_labels[mask]
        y_true = cond_labs == cond1
        n1 = y_true.sum()
        n2 = mask.sum() - n1

    if n1 > 0 and n2 > 0:
        for i in range(shape[1]):
            if cond_labels is None:
                exprs[:] = 0.0
                exprs[
                    csc_mat.indices[csc_mat.indptr[i] : csc_mat.indptr[i + 1]]
                ] = csc_mat.data[csc_mat.indptr[i] : csc_mat.indptr[i + 1]]
            else:
                exprs = csc_mat[mask, i].toarray()[:, 0]

            fpr, tpr, thresholds = sm.roc_curve(y_true, exprs)
            auroc[i] = sm.auc(fpr, tpr)

            # precision, recall, thresholds = sm.precision_recall_curve(y_true, exprs)
            # aupr[i] = sm.auc(recall, precision)

    df = pd.DataFrame(
        {
            "auroc:{0}".format(clust_id): auroc,
            # "aupr:{0}".format(clust_id): aupr,
        },
        index=gene_names,
    )

    if verbose:
        logger.info("calc_auc finished for cluster {0}.".format(clust_id))

    return df


def calculate_auc_values(
    Xc: csc_matrix,
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    n_jobs: int,
    temp_folder: str,
    verbose: bool,
) -> List[pd.DataFrame]:
    """ Calculate AUROC values, triggering calc_auc in parallel
    """
    start = time.perf_counter()

    result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
        delayed(calc_auc)(
            clust_id,
            Xc.data,
            Xc.indices,
            Xc.indptr,
            Xc.shape,
            cluster_labels,
            cond_labels,
            gene_names,
            verbose,
        )
        for clust_id in cluster_labels.categories
    )

    end = time.perf_counter()
    if verbose:
        logger.info(
            "AUROC values are calculated. Time spent = {:.2f}s.".format(end - start)
        )

    return result_list


def calc_t(
    clust_id: str,
    data: List[float],
    indices: List[int],
    indptr: List[int],
    shape: Tuple[int, int],
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    sum_vec: List[float],
    sum2_vec: List[float],
    verbose: bool,
) -> pd.DataFrame:
    """ Calcualte Welch's t-test for one cluster
    """

    # recover sparse matrix
    mat = csr_matrix((data, indices, indptr), shape=shape)
    pvals = np.full(shape[1], 1.0)
    tscores = np.full(shape[1], 0)
    mask = cluster_labels == clust_id
    mat_clust = mat[mask]

    if cond_labels is None:
        n1 = mat_clust.shape[0]
        n2 = shape[0] - n1

        if n1 > 1 and n2 > 1:
            sum_clu = mat_clust.sum(axis=0).A1
            mean1 = sum_clu / n1
            mean2 = (sum_vec - sum_clu) / n2
            mean2[mean2 < 0.0] = 0.0

            sum2_clu = mat_clust.power(2).sum(axis=0).A1
            s1sqr = (sum2_clu - n1 * (mean1 ** 2)) / (n1 - 1)
            s2sqr = ((sum2_vec - sum2_clu) - n2 * (mean2 ** 2)) / (n2 - 1)
            s2sqr[s2sqr < 0.0] = 0.0
    else:
        cond1 = cond_labels.categories[0]
        cond_labs = cond_labels[mask]
        mask2 = cond_labs == cond1

        mat_cond1 = mat_clust[mask2]
        mat_cond2 = mat_clust[~mask2]
        n1 = mat_cond1.shape[0]
        n2 = mat_cond2.shape[0]

        if n1 > 1 and n2 > 1:
            mean1 = mat_cond1.mean(axis=0).A1
            psum1 = mat_cond1.power(2).sum(axis=0).A1
            s1sqr = (psum1 - n1 * (mean1 ** 2)) / (n1 - 1)

            mean2 = mat_cond2.mean(axis=0).A1
            psum2 = mat_cond2.power(2).sum(axis=0).A1
            s2sqr = (psum2 - n2 * (mean2 ** 2)) / (n2 - 1)

    if n1 > 1 and n2 > 1:
        import scipy.stats as ss

        var_est = s1sqr / n1 + s2sqr / n2
        idx = var_est > 0.0
        if idx.sum() > 0:
            tscore = (mean1[idx] - mean2[idx]) / np.sqrt(var_est[idx])
            v = (var_est[idx] ** 2) / (
                (s1sqr[idx] / n1) ** 2 / (n1 - 1) + (s2sqr[idx] / n2) ** 2 / (n2 - 1)
            )
            pvals[idx] = ss.t.sf(np.fabs(tscore), v) * 2.0  # two-sided
            tscores[idx] = tscore
    passed, qvals = fdr(pvals)

    df = pd.DataFrame(
        {
            "t_pval:{0}".format(clust_id): pvals.astype(np.float32),
            "t_qval:{0}".format(clust_id): qvals.astype(np.float32),
            "t_score:{0}".format(clust_id): tscores.astype(np.float32),
        },
        index=gene_names,
    )

    if verbose:
        logger.info("calc_t finished for cluster {0}.".format(clust_id))

    return df


def t_test(
    X: csr_matrix,
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    n_jobs: int,
    temp_folder: str,
    verbose: bool,
) -> List[pd.DataFrame]:
    """ Run Welch's t-test, triggering calc_t in parallel
    """
    start = time.perf_counter()

    sum_vec = sum2_vec = None
    if cond_labels is None:
        sum_vec = X.sum(axis=0).A1
        sum2_vec = X.power(2).sum(axis=0).A1

    result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
        delayed(calc_t)(
            clust_id,
            X.data,
            X.indices,
            X.indptr,
            X.shape,
            cluster_labels,
            cond_labels,
            gene_names,
            sum_vec,
            sum2_vec,
            verbose,
        )
        for clust_id in cluster_labels.categories
    )

    end = time.perf_counter()
    if verbose:
        logger.info("Welch's t-test is done. Time spent = {:.2f}s.".format(end - start))

    return result_list


def calc_fisher(
    clust_id: str,
    data: List[float],
    indices: List[int],
    indptr: List[int],
    shape: Tuple[int, int],
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    cnt_vec: List[int],
    verbose: bool,
) -> pd.DataFrame:
    """ Calcualte Fisher's exact test for one cluster
    """
    import fisher

    # recover sparse matrix
    mat = csr_matrix((data, indices, indptr), shape=shape)
    mask = cluster_labels == clust_id
    mat_clust = mat[mask]

    if cond_labels is None:
        n1 = mat_clust.shape[0]
        n2 = shape[0] - n1

        a_true = mat_clust.getnnz(axis=0).astype(np.uint)
        a_false = n1 - a_true
        b_true = cnt_vec.astype(np.uint) - a_true
        b_false = n2 - b_true
    else:
        cond1 = cond_labels.categories[0]
        cond_labs = cond_labels[mask]
        mask2 = cond_labs == cond1

        mat_cond1 = mat_clust[mask2]
        mat_cond2 = mat_clust[~mask2]
        n1 = mat_cond1.shape[0]
        n2 = mat_cond2.shape[0]

        a_true = mat_cond1.getnnz(axis=0).astype(np.uint)
        a_false = n1 - a_true
        b_true = mat_cond2.getnnz(axis=0).astype(np.uint)
        b_false = n2 - b_true

    pvals = fisher.pvalue_npy(a_true, a_false, b_true, b_false)[2]
    passed, qvals = fdr(pvals)

    df = pd.DataFrame(
        {
            "fisher_pval:{0}".format(clust_id): pvals.astype(np.float32),
            "fisher_qval:{0}".format(clust_id): qvals.astype(np.float32),
        },
        index=gene_names,
    )

    if verbose:
        logger.info("calc_fisher finished for cluster {0}.".format(clust_id))

    return df


def fisher_test(
    X: csr_matrix,
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    n_jobs: int,
    temp_folder: str,
    verbose: bool,
) -> List[pd.DataFrame]:
    """ Run Fisher's exact test, triggering calc_fisher in parallel
    """
    start = time.perf_counter()

    cnt_vec = None
    if cond_labels is None:
        cnt_vec = X.getnnz(axis=0)

    result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
        delayed(calc_fisher)(
            clust_id,
            X.data,
            X.indices,
            X.indptr,
            X.shape,
            cluster_labels,
            cond_labels,
            gene_names,
            cnt_vec,
            verbose,
        )
        for clust_id in cluster_labels.categories
    )

    end = time.perf_counter()
    if verbose:
        logger.info(
            "Fisher's exact test is done. Time spent = {:.2f}s.".format(end - start)
        )

    return result_list


def calc_mwu(
    clust_id: str,
    data: List[float],
    indices: List[int],
    indptr: List[int],
    shape: Tuple[int, int],
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    verbose: bool,
) -> pd.DataFrame:
    """ Run Mann-Whitney U test for one cluster
    """

    csc_mat = csc_matrix((data, indices, indptr), shape=shape)
    U_stats = np.zeros(shape[1], dtype=np.float32)
    pvals = np.full(shape[1], 1.0)
    mask = cluster_labels == clust_id

    if cond_labels is None:
        exprs = np.zeros(shape[0])
        idx_x = mask
        idx_y = ~idx_x
    else:
        exprs = None
        cond1 = cond_labels.categories[0]
        cond_labs = cond_labels[mask]
        idx_x = cond_labs == cond1
        idx_y = ~idx_x

    n1 = idx_x.sum()
    n2 = idx_y.sum()

    if n1 > 0 and n2 > 0:
        import scipy.stats as ss

        for i in range(shape[1]):
            if cond_labels is None:
                if csc_mat.indptr[i + 1] - csc_mat.indptr[i] > 0:
                    exprs[:] = 0.0
                    exprs[
                        csc_mat.indices[csc_mat.indptr[i] : csc_mat.indptr[i + 1]]
                    ] = csc_mat.data[csc_mat.indptr[i] : csc_mat.indptr[i + 1]]
                    U_stats[i], pvals[i] = ss.mannwhitneyu(
                        exprs[idx_x], exprs[idx_y], alternative="two-sided"
                    )
            else:
                tmp_mat = csc_mat[mask, i]
                if tmp_mat.data.size > 0:
                    exprs = tmp_mat.toarray()[:, 0]
                    U_stats[i], pvals[i] = ss.mannwhitneyu(
                        exprs[idx_x], exprs[idx_y], alternative="two-sided"
                    )

    passed, qvals = fdr(pvals)

    df = pd.DataFrame(
        {
            "mwu_U:{0}".format(clust_id): U_stats.astype(np.float32),
            "mwu_pval:{0}".format(clust_id): pvals.astype(np.float32),
            "mwu_qval:{0}".format(clust_id): qvals.astype(np.float32),
        },
        index=gene_names,
    )

    if verbose:
        logger.info("calc_mwu finished for cluster {0}.".format(clust_id))

    return df


def mwu_test(
    Xc: csc_matrix,
    cluster_labels: List[str],
    cond_labels: List[str],
    gene_names: List[str],
    n_jobs: int,
    temp_folder: str,
    verbose: bool,
) -> List[pd.DataFrame]:
    """ Run Mann-Whitney U test, triggering calc_mwu in parallel
    """
    start = time.perf_counter()

    result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
        delayed(calc_mwu)(
            clust_id,
            Xc.data,
            Xc.indices,
            Xc.indptr,
            Xc.shape,
            cluster_labels,
            cond_labels,
            gene_names,
            verbose,
        )
        for clust_id in cluster_labels.categories
    )

    end = time.perf_counter()
    if verbose:
        logger.info(
            "Mann-Whitney U test is done. Time spent = {:.2f}s.".format(end - start)
        )

    return result_list


def organize_results(results: List[List[pd.DataFrame]]) -> pd.DataFrame:
    """ Concatenate resulting dataframes into one big dataframe
    """
    m = len(results)
    n = len(results[0])
    reslist = []

    for i in range(n):
        for j in range(m):
            reslist.append(results[j][i])

    df = pd.concat(reslist, axis=1)
    return df


def de_analysis(
    data: Union[MultimodalData, UnimodalData, AnnData],
    cluster: str,
    condition: str = None,
    subset: str = None,
    result_key: str = "de_res",
    n_jobs: int = -1,
    auc: bool = True,
    t: bool = True,
    fisher: bool = False,
    mwu: bool = False,
    temp_folder: str = None,
    verbose: bool = True,
) -> None:
    """Perform Differential Expression (DE) Analysis on data.

    Parameters
    ----------
    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``
        Data matrix with rows for cells and columns for genes.

    cluster: ``str``
        Cluster labels used in DE analysis. Must exist in ``data.obs``.

    condition: ``str``, optional, default: ``None``
        Sample attribute used as condition in DE analysis. If ``None``, no condition is considered; otherwise, must exist in ``data.obs``.

    subset: ``str``, optional, default: ``None``
        Perform DE analysis on only a subset of cluster IDs. Cluster ID subset is specified as ``"clust_id1,clust_id2,...,clust_idn"``, where all IDs must exist in ``data.obs[cluster]``.

    result_key: ``str``, optional, default: ``"de_res"``
        Key name of DE analysis result stored.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    auc: ``bool``, optional, default: ``True``
        If ``True``, calculate area under ROC (AUROC) and area under Precision-Recall (AUPR).

    t: ``bool``, optional, default: ``True``
        If ``True``, calculate Welch's t test.

    fisher: ``bool``, optional, default: ``False``
        If ``True``, calculate Fisher's exact test.

    mwu: ``bool``, optional, default: ``False``
        If ``True``, calculate Mann-Whitney U test.

    temp_folder: ``str``, optional, default: ``None``
        Joblib temporary folder for memmapping numpy arrays.

    verbose: ``bool``, optional, default: ``True``
        If ``True``, show detailed intermediate output.

    Returns
    -------
    ``None``

    Update ``data.varm``:
        ``data.varm[result_key]``: DE analysis result.

    Examples
    --------
    >>> pg.de_analysis(data, cluster = 'spectral_leiden_labels')

    subset: a comma-separated list of cluster labels. Then de will be performed only on these subsets.
    """
    start = time.perf_counter()

    if cluster not in data.obs:
        raise ValueError("Cannot find cluster label!")
    cluster_labels = data.obs[cluster].values
    if not isinstance(cluster_labels, pd.Categorical):
        cluster_labels = pd.Categorical(cluster_labels)

    cond_labels = None
    if condition is not None:
        if condition not in data.obs:
            raise ValueError("Cannot find condition!")
        cond_labels = data.obs[condition].values
        if not isinstance(cond_labels, pd.Categorical):
            cond_labels = pd.Categorical(cond_labels)
        if cond_labels.categories.size != 2:
            raise ValueError(
                "Number of distinct values in Condition is not equal to 2!"
            )

    X = data.X if isinstance(data.X, csr_matrix) else data.X[:]
    X.eliminate_zeros()  # In case there is any extra zeros

    if subset is not None:
        # subset data for de analysis
        subset = np.array(subset.split(","))
        idx_s = np.isin(subset, cluster_labels.categories.values)
        if idx_s.sum() < subset.size:
            raise ValueError(
                "These cluster labels do not exist: " + ",".join(subset[~idx_s]) + "!"
            )

        idx = np.isin(cluster_labels, subset)
        cluster_labels = cluster_labels[idx]
        if cond_labels is not None:
            cond_labels = cond_labels[idx]
        X = X[idx]

    n_jobs = effective_n_jobs(n_jobs)
    gene_names = data.var_names

    results = []
    results.append(
        collect_basic_statistics(
            X, cluster_labels, cond_labels, gene_names, n_jobs, temp_folder, verbose
        )
    )

    Xc = None
    if auc or mwu:
        t1 = time.perf_counter()
        Xc = X.tocsc()
        if verbose:
            logger.info(
                "Converting X to csc_matrix is done. Time spent = {:.2f}s.".format(
                    time.perf_counter() - t1
                )
            )

    if auc:
        results.append(
            calculate_auc_values(
                Xc,
                cluster_labels,
                cond_labels,
                gene_names,
                n_jobs,
                temp_folder,
                verbose,
            )
        )

    if t:
        results.append(
            t_test(
                X, cluster_labels, cond_labels, gene_names, n_jobs, temp_folder, verbose
            )
        )

    if fisher:
        results.append(
            fisher_test(
                X, cluster_labels, cond_labels, gene_names, n_jobs, temp_folder, verbose
            )
        )

    if mwu:
        results.append(
            mwu_test(
                Xc,
                cluster_labels,
                cond_labels,
                gene_names,
                n_jobs,
                temp_folder,
                verbose,
            )
        )

    df = organize_results(results)
    data.varm[result_key] = df.to_records(index=False)

    end = time.perf_counter()
    logger.info(
        "Differential expression analysis is finished. Time spent = {:.2f}s.".format(
            end - start
        )
    )


def get_valid_gene_index(n: int, df: pd.DataFrame, alpha: float) -> List[bool]:
    """ get genes that are DE for at least one test. If no DE tests, all genes are valid.
    """
    idx = np.zeros(n, dtype=bool)
    has_test = False
    for qval in ["t_qval", "fisher_qval", "mwu_qval"]:
        if qval in df.columns:
            idx = idx | (df[qval].values <= alpha)
            has_test = True

    if not has_test:
        idx = np.ones(n, dtype=bool)

    return idx


def get_sort_key(sort_by: List[str], col_names: List[str], direction: str):
    """ direction is either up or down; if down, do not use aupr
    """
    for key in sort_by:
        if key in col_names:
            return key
    raise ValueError("No valid key!")


def markers(
    data: Union[MultimodalData, UnimodalData, AnnData],
    head: int = None,
    de_key: str = "de_res",
    sort_by: str = "auroc,WAD_score",
    alpha: float = 0.05,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """

    Parameters
    ----------
    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``
        Data matrix with rows for cells and columns for genes.

    head: ``int``, optional, default: ``None``
        List only top ``head`` genes for each cluster. If ``None``, show any DE genes.

    de_key: ``str``, optional, default, ``de_res``
        Keyword of DE result stored in ``data.varm``.

    sort_by: ``str``, optional, default: ``"auroc,WAD_score"``
        Sort the resulting marker dictionary by ``auroc`` and ``WAD_score``.

    alpha: ``float``, optional, default: ``0.05``
        q-value threshold for getting valid DE genes. Only those with q-value of any test below ``alpha`` are significant, and thus considered as DE genes.

    Returns
    -------
    results: ``Dict[str, Dict[str, pd.DataFrame]]``
        A Python dictionary containing markers in structure ``dict[cluster_id]['up' or 'down'][dataframe]``.

    Examples
    --------
    >>> marker_dict = pg.markers(data)
    """
    if de_key not in data.varm.keys():
        raise ValueError("Please run de_analysis first!")

    sort_by = sort_by.split(",")

    clust2cols = defaultdict(list)
    rec_array = data.varm[de_key]
    for name in rec_array.dtype.names:
        col_name, sep, clust_id = name.partition(":")
        clust2cols[clust_id].append(col_name)

    results = defaultdict(dict)
    for clust_id, col_names in clust2cols.items():
        rec_names = [x + ":" + clust_id for x in col_names]
        df = pd.DataFrame(data=rec_array[rec_names], index=data.var_names.copy())
        df.columns = col_names
        df.index.name = "feature"

        idx = get_valid_gene_index(data.shape[1], df, alpha)

        idx_up = idx & (df["log_fold_change"].values > 0)
        df_up = df.loc[idx_up].sort_values(
            by=get_sort_key(sort_by, col_names, "up"), ascending=False, inplace=False
        )
        results[clust_id]["up"] = pd.DataFrame(
            df_up if head is None else df_up.iloc[0:head]
        )

        idx_down = idx & (df["log_fold_change"].values < 0)
        df_down = df.loc[idx_down].sort_values(
            by=get_sort_key(sort_by, col_names, "down"), ascending=True, inplace=False
        )
        results[clust_id]["down"] = pd.DataFrame(
            df_down if head is None else df_down.iloc[0:head]
        )

    return results


def write_results_to_excel(
    results: Dict[str, Dict[str, pd.DataFrame]], output_file: str, ndigits: int = 3
) -> None:
    """ Write results into Excel workbook.

    Parameters
    ----------
    results: ``Dict[str, Dict[str, pd.DataFrame]]``
        DE marker dictionary generated by ``pg.markers``.

    output_file: ``str``
        File name to which the marker dictionary is written.

    ndigits: ``int``, optional, default: ``3``
        Round non p-values and q-values to ``ndigits`` after decimal point in the excel.

    Returns
    -------
    ``None``

    Marker information is written to file with name ``output_file``.

    Examples
    --------
    >>> pg.write_results_to_excel(marker_dict, "result.de.xlsx")
    """
    start = time.perf_counter()

    import xlsxwriter
    from natsort import natsorted

    def format_short_output_cols(
        df_orig: pd.DataFrame, ndigits: int = 3
    ) -> pd.DataFrame:
        """ Round related float columns to ndigits decimal points.
        """
        df = pd.DataFrame(df_orig)

        cols = []
        for name in df.columns:
            if (not name.endswith("pval")) and (not name.endswith("qval")):
                cols.append(name)

        df.loc[:, cols] = df.loc[:, cols].round(ndigits)
        return df

    def add_worksheet(
        workbook: "workbook", df_orig: pd.DataFrame, sheet_name: str
    ) -> None:
        """ Add one worksheet with content as df
        """
        df = format_short_output_cols(df_orig)
        df.reset_index(inplace=True)
        worksheet = workbook.add_worksheet(name=sheet_name)

        if df.shape[0] > 0:
            worksheet.add_table(
                0,
                0,
                df.index.size,
                df.columns.size - 1,
                {
                    "data": df.to_numpy(),
                    "style": "Table Style Light 1",
                    "first_column": True,
                    "header_row": True,
                    "columns": [{"header": x} for x in df.columns.values],
                },
            )
        else:
            worksheet.write_row(0, 0, df.columns.values)

    workbook = xlsxwriter.Workbook(output_file, {"nan_inf_to_errors": True})
    workbook.formats[0].set_font_size(9)

    for clust_id in natsorted(results.keys()):
        add_worksheet(workbook, results[clust_id]["up"], "up " + clust_id)
        add_worksheet(workbook, results[clust_id]["down"], "down " + clust_id)

    workbook.close()

    end = time.perf_counter()
    logger.info(
        "Excel spreadsheet is written. Time spent = {:.2f}s.".format(end - start)
    )


@timer(logger=logger)
def run_de_analysis(
    input_file: str,
    output_excel_file: str,
    cluster: str,
    result_key: str = "de_res",
    n_jobs: int = -1,
    auc: bool = True,
    t: bool = True,
    fisher: bool = False,
    mwu: bool = False,
    temp_folder: str = None,
    verbose: bool = True,
    alpha: float = 0.05,
    ndigits: int = 3,
) -> None:
    """ For command line only
    """

    from pegasusio import read_input, write_output

    data = read_input(input_file, mode='r')

    de_analysis(
        data,
        cluster,
        result_key=result_key,
        n_jobs=n_jobs,
        auc=auc,
        t=t,
        fisher=fisher,
        mwu=mwu,
        temp_folder=temp_folder,
        verbose=verbose,
    )

    write_output(data, input_file)
    logger.info(
        "Differential expression results are written to varm/{} in zarr file.".format(
            result_key
        )
    )

    results = markers(data, de_key=result_key, alpha=alpha)

    write_results_to_excel(results, output_excel_file, ndigits=ndigits)
