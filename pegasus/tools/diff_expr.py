import time
import numpy as np
import pandas as pd

from pandas.api.types import is_categorical_dtype
from scipy.sparse import csr_matrix
from statsmodels.stats.multitest import fdrcorrection as fdr
from joblib import Parallel, delayed, parallel_backend

from typing import List, Tuple, Dict, Union, Optional

import logging
logger = logging.getLogger(__name__)

from anndata import AnnData
from pegasusio import timer, MultimodalData, UnimodalData
from pegasus.tools import eff_n_jobs



def _calc_qvals(
    nclust: int,
    pvals: np.ndarray,
    first_j: int,
    second_j: int,
) -> np.ndarray:
    """ Calculate FDR
    """
    qvals = np.zeros(pvals.shape, dtype = np.float32)
    if second_j > 0:
        _, qval = fdr(pvals[:, first_j])
        qvals[:, first_j] = qvals[:, second_j] = qval
    else:
        for j in range(nclust):
            _, qvals[:, j] = fdr(pvals[:, j])
    return qvals


def _de_test(
    X: csr_matrix,
    cluster_labels: pd.Categorical,
    gene_names: List[str],
    n_jobs: int,
    t: Optional[bool] = False,
    fisher: Optional[bool] = False,
    temp_folder: Optional[str] = None,
    verbose: Optional[bool] = True,
) -> pd.DataFrame:
    """ Collect sufficient statistics, run Mann-Whitney U test, calculate auroc (triggering diff_expr_utils.calc_mwu in parallel), optionally run Welch's T test and Fisher's Exact test (in parallel).
    """
    from pegasus.cylib.de_utils import csr_to_csc, calc_mwu, calc_stat

    start = time.perf_counter()

    ords = np.argsort(cluster_labels.codes)
    data, indices, indptr = csr_to_csc(X.data, X.indices, X.indptr, X.shape[0], X.shape[1], ords)
    cluster_cnts = cluster_labels.value_counts()
    n1arr = cluster_cnts.values
    n2arr = X.shape[0] - n1arr
    cluster_cumsum = cluster_cnts.cumsum().values
    nclust = n1arr.size

    first_j = second_j = -1
    posvec = np.where(n1arr > 0)[0]
    if len(posvec) == 2:
        first_j = posvec[0]
        second_j = posvec[1]


    if verbose:
        end = time.perf_counter()
        logger.info(f"CSR matrix is converted to CSC matrix. Time spent = {end - start:.4f}s.")
        start = end
        # logger.info(f"Preparation (including converting X to csc_matrix format) for MWU test is finished. Time spent = {time.perf_counter() - start:.2f}s.")


    ngene = X.shape[1]
    quotient = ngene // n_jobs
    residue = ngene % n_jobs
    intervals = []
    start_pos = end_pos = 0
    for i in range(n_jobs):
        end_pos = start_pos + quotient + (i < residue)
        if end_pos == start_pos:
            break
        intervals.append((start_pos, end_pos))
        start_pos = end_pos

    with parallel_backend("loky", inner_max_num_threads=1):
        result_list = Parallel(n_jobs=len(intervals), temp_folder=temp_folder)(
            delayed(calc_mwu)(
                start_pos,
                end_pos,
                data,
                indices,
                indptr,
                n1arr,
                n2arr,
                cluster_cumsum,
                first_j,
                second_j,
                verbose,
            )
            for start_pos, end_pos in intervals
        )

    Ulist = []
    plist = []
    alist = []
    for U_stats, pvals, aurocs in result_list:
        Ulist.append(U_stats)
        plist.append(pvals)
        alist.append(aurocs)

    U_stats = np.concatenate(Ulist, axis = 0)
    pvals = np.concatenate(plist, axis = 0)
    aurocs = np.concatenate(alist, axis = 0)
    qvals = _calc_qvals(nclust, pvals, first_j, second_j)

    dfU = pd.DataFrame(U_stats, index = gene_names, columns = [f"{x}:mwu_U" for x in cluster_labels.categories])
    dfUp = pd.DataFrame(pvals, index = gene_names, columns = [f"{x}:mwu_pval" for x in cluster_labels.categories])
    dfUq = pd.DataFrame(qvals, index = gene_names, columns = [f"{x}:mwu_qval" for x in cluster_labels.categories])
    dfUa = pd.DataFrame(aurocs, index = gene_names, columns = [f"{x}:auroc" for x in cluster_labels.categories])

    if verbose:
        end = time.perf_counter()
        logger.info(f"MWU test and AUROC calculation are finished. Time spent = {end - start:.4f}s.")
        start = end

    # basic statistics and optional t test and fisher test
    results = calc_stat(data, indices, indptr, n1arr, n2arr, cluster_cumsum, first_j, second_j, t, fisher, verbose)
    dfl2M = pd.DataFrame(results[0][0], index = gene_names, columns = [f"{x}:log2Mean" for x in cluster_labels.categories])
    dfl2Mo = pd.DataFrame(results[0][1], index = gene_names, columns = [f"{x}:log2Mean_other" for x in cluster_labels.categories])
    dfl2FC = pd.DataFrame(results[0][2], index = gene_names, columns = [f"{x}:log2FC" for x in cluster_labels.categories])
    dfpct = pd.DataFrame(results[0][3], index = gene_names, columns = [f"{x}:percentage" for x in cluster_labels.categories])
    dfpcto = pd.DataFrame(results[0][4], index = gene_names, columns = [f"{x}:percentage_other" for x in cluster_labels.categories])
    dfpfc = pd.DataFrame(results[0][5], index = gene_names, columns = [f"{x}:percentage_fold_change" for x in cluster_labels.categories])

    df_list = [dfl2M, dfl2Mo, dfl2FC, dfpct, dfpcto, dfpfc, dfUa, dfU, dfUp, dfUq]

    if verbose:
        end = time.perf_counter()
        logger.info(f"Sufficient statistics are collected. Time spent = {end - start:.4f}s.")
        start = end

    if t:
        qvals = _calc_qvals(nclust, results[1][1], first_j, second_j)
        dft = pd.DataFrame(results[1][0], index = gene_names, columns = [f"{x}:t_tstat" for x in cluster_labels.categories])
        dftp = pd.DataFrame(results[1][1], index = gene_names, columns = [f"{x}:t_pval" for x in cluster_labels.categories])
        dftq = pd.DataFrame(qvals, index = gene_names, columns = [f"{x}:t_qval" for x in cluster_labels.categories])
        df_list.extend([dft, dftp, dftq])

        if verbose:
            end = time.perf_counter()
            logger.info(f"Welch's t-test is finished. Time spent = {end - start:.4f}s.")
            start = end

    if fisher:
        from pegasus.cylib.cfisher import fisher_exact
        a_true, a_false, b_true, b_false = results[1 if not t else 2]

        oddsratios = np.zeros((ngene, n1arr.size), dtype = np.float32)
        pvals = np.ones((ngene, n1arr.size), dtype = np.float32)

        if second_j > 0:
            oddsratio, pval = fisher_exact(a_true[first_j], a_false[first_j], b_true[first_j], b_false[first_j])
            oddsratios[:, first_j] = oddsratio
            idx1 = oddsratio > 0.0
            idx2 = oddsratio < 1e30
            oddsratios[idx1 & idx2, second_j] = 1.0 / oddsratio[idx1 & idx2]
            oddsratios[~idx1] = 1e30
            pvals[:, first_j] = pvals[:, second_j] = pval
        else:
            with parallel_backend("loky", inner_max_num_threads=1):
                result_list = Parallel(n_jobs=n_jobs, temp_folder=temp_folder)(
                    delayed(fisher_exact)(
                        a_true[i],
                        a_false[i],
                        b_true[i],
                        b_false[i],
                    )
                    for i in posvec
                )

            for i in range(posvec.size):
                oddsratios[:, posvec[i]] = result_list[i][0]
                pvals[:, posvec[i]] = result_list[i][1]

        qvals = _calc_qvals(nclust, pvals, first_j, second_j)
        dff = pd.DataFrame(oddsratios, index = gene_names, columns = [f"{x}:fisher_oddsratio" for x in cluster_labels.categories])
        dffp = pd.DataFrame(pvals, index = gene_names, columns = [f"{x}:fisher_pval" for x in cluster_labels.categories])
        dffq = pd.DataFrame(qvals, index = gene_names, columns = [f"{x}:fisher_qval" for x in cluster_labels.categories])
        df_list.extend([dff, dffp, dffq])

        if verbose:
            end = time.perf_counter()
            logger.info(f"Fisher's exact test is finished. Time spent = {end - start:.4f}s.")

    df = pd.concat(df_list, axis = 1)

    return df



def _perform_de_cond(
    clust_ids: List[str],
    cond_labels: pd.Categorical,
    gene_names: List[str],
    cond_n1arr_list: List[List[int]],
    cond_n2arr_list: List[List[int]],
    cond_cumsum_list: List[List[int]],
    data_list: List[List[float]],
    indices_list: List[List[int]],
    indptr_list: List[List[int]],
    t: bool,
    fisher: bool,
    verbose: bool,
) -> List[pd.DataFrame]:
    """ Run DE test for clusters specified. In each cluster, perform one vs. rest for the condition
    """
    df_res_list = []

    from pegasus.cylib.de_utils import calc_mwu, calc_stat

    ngene = indptr_list[0].size - 1
    for i, clust_id in enumerate(clust_ids):
        nclust = cond_n1arr_list[i].size

        first_j = second_j = -1
        posvec = np.where(cond_n1arr_list[i] > 0)[0]
        if len(posvec) == 2:
            first_j = posvec[0]
            second_j = posvec[1]

        U_stats, pvals, aurocs = calc_mwu(0, ngene, data_list[i], indices_list[i], indptr_list[i], cond_n1arr_list[i], cond_n2arr_list[i], cond_cumsum_list[i], first_j, second_j, False)
        qvals = _calc_qvals(nclust, pvals, first_j, second_j)

        dfU = pd.DataFrame(U_stats, index = gene_names, columns = [f"{clust_id}:{x}:mwu_U" for x in cond_labels.categories])
        dfUp = pd.DataFrame(pvals, index = gene_names, columns = [f"{clust_id}:{x}:mwu_pval" for x in cond_labels.categories])
        dfUq = pd.DataFrame(qvals, index = gene_names, columns = [f"{clust_id}:{x}:mwu_qval" for x in cond_labels.categories])
        dfUa = pd.DataFrame(aurocs, index = gene_names, columns = [f"{clust_id}:{x}:auroc" for x in cond_labels.categories])

        results = calc_stat(data_list[i], indices_list[i], indptr_list[i], cond_n1arr_list[i], cond_n2arr_list[i], cond_cumsum_list[i], first_j, second_j, t, fisher, False)
        dfl2M = pd.DataFrame(results[0][0], index = gene_names, columns = [f"{clust_id}:{x}:log2Mean" for x in cond_labels.categories])
        dfl2Mo = pd.DataFrame(results[0][1], index = gene_names, columns = [f"{clust_id}:{x}:log2Mean_other" for x in cond_labels.categories])
        dfl2FC = pd.DataFrame(results[0][2], index = gene_names, columns = [f"{clust_id}:{x}:log2FC" for x in cond_labels.categories])
        dfpct = pd.DataFrame(results[0][3], index = gene_names, columns = [f"{clust_id}:{x}:percentage" for x in cond_labels.categories])
        dfpcto = pd.DataFrame(results[0][4], index = gene_names, columns = [f"{clust_id}:{x}:percentage_other" for x in cond_labels.categories])
        dfpfc = pd.DataFrame(results[0][5], index = gene_names, columns = [f"{clust_id}:{x}:percentage_fold_change" for x in cond_labels.categories])

        df_list = [dfl2M, dfl2Mo, dfl2FC, dfpct, dfpcto, dfpfc, dfUa, dfU, dfUp, dfUq]

        if t:
            qvals = _calc_qvals(nclust, results[1][1], first_j, second_j)
            dft = pd.DataFrame(results[1][0], index = gene_names, columns = [f"{clust_id}:{x}:t_tstat" for x in cond_labels.categories])
            dftp = pd.DataFrame(results[1][1], index = gene_names, columns = [f"{clust_id}:{x}:t_pval" for x in cond_labels.categories])
            dftq = pd.DataFrame(qvals, index = gene_names, columns = [f"{clust_id}:{x}:t_qval" for x in cond_labels.categories])
            df_list.extend([dft, dftp, dftq])

        if fisher:
            a_true, a_false, b_true, b_false = results[1 if not t else 2]

            from pegasus.cylib.cfisher import fisher_exact

            oddsratios = np.zeros((ngene, nclust), dtype = np.float32)
            pvals = np.ones((ngene, nclust), dtype = np.float32)

            if second_j > 0:
                oddsratio, pval = fisher_exact(a_true[first_j], a_false[first_j], b_true[first_j], b_false[first_j])
                oddsratios[:, first_j] = oddsratio
                idx1 = oddsratio > 0.0
                idx2 = oddsratio < 1e30
                oddsratios[idx1 & idx2, second_j] = 1.0 / oddsratio[idx1 & idx2]
                oddsratios[~idx1] = 1e30
                pvals[:, first_j] = pvals[:, second_j] = pval
            else:
                for j in posvec:
                    oddsratios[:, j], pvals[:, j] = fisher_exact(a_true[j], a_false[j], b_true[j], b_false[j])

            qvals = _calc_qvals(nclust, pvals, first_j, second_j)
            dff = pd.DataFrame(oddsratios, index = gene_names, columns = [f"{clust_id}:{x}:fisher_oddsratio" for x in cond_labels.categories])
            dffp = pd.DataFrame(pvals, index = gene_names, columns = [f"{clust_id}:{x}:fisher_pval" for x in cond_labels.categories])
            dffq = pd.DataFrame(qvals, index = gene_names, columns = [f"{clust_id}:{x}:fisher_qval" for x in cond_labels.categories])
            df_list.extend([dff, dffp, dffq])

        df_res_list.append(pd.concat(df_list, axis = 1))

    if verbose:
        clust_ids_str = ','.join([str(x) for x in clust_ids])
        print(f"_perform_de_cond finished for clusters {clust_ids_str}.")

    return df_res_list


def _de_test_cond(
    X: csr_matrix,
    cluster_labels: pd.Categorical,
    cond_labels: pd.Categorical,
    gene_names: List[str],
    n_jobs: int,
    t: Optional[bool] = False,
    fisher: Optional[bool] = False,
    temp_folder: Optional[str] = None,
    verbose: Optional[bool] = True,
) -> List[pd.DataFrame]:
    """ Collect sufficient statistics, run Mann-Whitney U test, calculate auroc, optionally run Welch's T test and Fisher's Exact test on each cluster (DE analysis is based on cond_labels)
    """
    start = time.perf_counter()

    clust_cond = np.array(list(zip(cluster_labels.codes, cond_labels.codes)), dtype = [("clust", "<i4"), ("cond", "<i4")])
    ords = np.argsort(clust_cond)
    df_cross = pd.crosstab(cluster_labels.codes, cond_labels.codes)
    cluster_cnts = df_cross.values.sum(axis = 1)
    cluster_cumsum = cluster_cnts.cumsum()

    from pegasus.cylib.de_utils import csr_to_csc_cond

    data_list, indices_list, indptrs = csr_to_csc_cond(X.data, X.indices, X.indptr, X.shape[1], ords, cluster_cumsum)

    if verbose:
        end = time.perf_counter()
        logger.info(f"CSR matrix is converted to CSC matrix. Time spent = {end - start:.2f}s.")

    # assign from 0 to n and then n to 0; cluster ID is sorted in descending order with respect to number of cells
    ords = np.argsort(cluster_cnts)[::-1]
    nclust = cluster_cumsum.size
    neff = min(n_jobs, nclust)

    intervals = []
    datalists = []
    indiceslists = []
    for i in range(neff):
        intervals.append([])
        datalists.append([])
        indiceslists.append([])

    pos = 0
    sign = 1
    for i in range(nclust):
        intervals[pos].append(ords[i])
        datalists[pos].append(data_list[ords[i]])
        indiceslists[pos].append(indices_list[ords[i]])
        if pos + sign == neff:
            sign = -1
        elif pos + sign == -1:
            sign = 1
        else:
            pos += sign

    n1arr = df_cross.values
    n2arr = cluster_cnts.reshape(-1, 1) - n1arr
    cumsum = df_cross.cumsum(axis = 1).values

    with parallel_backend("loky", inner_max_num_threads=1):
        result_list = Parallel(n_jobs=neff, temp_folder=temp_folder)(
            delayed(_perform_de_cond)(
                cluster_labels.categories[intervals[i]],
                cond_labels,
                gene_names,
                n1arr[intervals[i]],
                n2arr[intervals[i]],
                cumsum[intervals[i]],
                datalists[i],
                indiceslists[i],
                indptrs[intervals[i]],
                t,
                fisher,
                verbose,
            )
            for i in range(neff)
        )

    df = pd.concat([x for y in result_list for x in y], axis = 1)

    return df


@timer(logger=logger)
def de_analysis(
    data: Union[MultimodalData, UnimodalData, AnnData],
    cluster: str,
    condition: Optional[str] = None,
    subset: Optional[List[str]] = None,
    de_key: Optional[str] = "de_res",
    n_jobs: Optional[int] = -1,
    t: Optional[bool] = False,
    fisher: Optional[bool] = False,
    temp_folder: Optional[str] = None,
    verbose: Optional[bool] = True,
) -> None:
    """Perform Differential Expression (DE) Analysis on data.

    The analysis considers one cluster at one time, comparing gene expression levels on cells
    within the cluster with all the others using a number of statistical tools, and determining
    up-regulated genes and down-regulated genes of the cluster.

    Mann-Whitney U test and AUROC are calculated by default. Welch's T test and Fisher's Exact test are optionally.

    The scalability performance on calculating all the test statistics is improved by the inspiration from `Presto <https://github.com/immunogenomics/presto>`_.

    Parameters
    ----------
    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``
        Data matrix with rows for cells and columns for genes.

    cluster: ``str``
        Cluster labels used in DE analysis. Must exist in ``data.obs``.

    condition: ``str``, optional, default: ``None``
        Sample attribute used as condition in DE analysis. If ``None``, no condition is considered; otherwise, must exist in ``data.obs``.
        If ``condition`` is used, the DE analysis will be performed on cells of each level of ``data.obs[condition]`` respectively, and collect the results after finishing.

    subset: ``List[str]``, optional, default: ``None``
        Perform DE analysis on only a subset of cluster IDs. Cluster ID subset is specified as a list of strings, such as ``[clust_1,clust_3,clust_5]``, where all IDs must exist in ``data.obs[cluster]``.

    de_key: ``str``, optional, default: ``"de_res"``
        Key name of DE analysis results stored.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all available threads.

    t: ``bool``, optional, default: ``True``
        If ``True``, calculate Welch's t test.

    fisher: ``bool``, optional, default: ``False``
        If ``True``, calculate Fisher's exact test.

    temp_folder: ``str``, optional, default: ``None``
        Joblib temporary folder for memmapping numpy arrays.

    verbose: ``bool``, optional, default: ``True``
        If ``True``, show detailed intermediate output.

    Returns
    -------
    ``None``

    Update ``data.varm``:
        ``data.varm[de_key]``: DE analysis result.

    Examples
    --------
    >>> pg.de_analysis(data, cluster='spectral_leiden_labels')
    >>> pg.de_analysis(data, cluster='louvain_labels', condition='anno')
    """
    if cluster not in data.obs:
        raise ValueError("Cannot find cluster label!")
    cluster_labels = data.obs[cluster].values
    if not is_categorical_dtype(cluster_labels):
        from natsort import natsorted
        cluster_labels = pd.Categorical(cluster_labels, natsorted(np.unique(cluster_labels)))

    cond_labels = None
    if condition is not None:
        if condition not in data.obs:
            raise ValueError("Cannot find condition!")
        cond_labels = data.obs[condition].values
        if not is_categorical_dtype(cond_labels):
            from natsort import natsorted
            cond_labels = pd.Categorical(cond_labels, natsorted(np.unique(cond_labels)))
        if cond_labels.categories.size < 2:
            raise ValueError("Number of conditions must be at least 2!")

    X = data.X if isinstance(data.X, csr_matrix) else csr_matrix(data.X) # If dense matrix, force it to be a csr_matrix

    if subset is not None:
        # subset data for de analysis
        subset = np.array(subset)
        idx_s = np.isin(subset, cluster_labels.categories.values)
        if idx_s.sum() < subset.size:
            raise ValueError(
                "These cluster labels do not exist: " + ",".join(subset[~idx_s]) + "!"
            )

        idx = np.isin(cluster_labels, subset)
        cluster_labels = pd.Categorical(cluster_labels[idx], categories = subset)
        if cond_labels is not None:
            cond_labels = cond_labels[idx]
        X = X[idx]

    if condition is not None:
        #Eliminate NaN rows from calculation
        idx_na = cond_labels.isna()
        if idx_na.sum() > 0:
            logger.warning("Detected NaN values in condition. Cells with NaN values are excluded from DE analysis.")
            idx_not_na = ~idx_na
            X = X[idx_not_na]
            cluster_labels = cluster_labels[idx_not_na]
            cond_labels = cond_labels[idx_not_na]

    n_jobs = eff_n_jobs(n_jobs)
    gene_names = data.var_names.values

    if cond_labels is None:
        df = _de_test(X, cluster_labels, gene_names, n_jobs, t, fisher, temp_folder, verbose)
    else:
        df = _de_test_cond(X, cluster_labels, cond_labels, gene_names, n_jobs, t, fisher, temp_folder, verbose)

    data.varm[de_key] = df.to_records(index=False)

    logger.info("Differential expression analysis is finished.")



def _assemble_df(res_dict: dict, rec_array: np.ndarray, prefix: str, col_names: List[str], gene_names: pd.Index, alpha: float, head: int):
    rec_names = [f"{prefix}:{x}" for x in col_names]
    df = pd.DataFrame(data=rec_array[rec_names], index=gene_names)
    df.columns = col_names

    idx = df["mwu_qval"] <= alpha
    idx_up = idx & (df["auroc"].values > 0.5)
    df_up = df.loc[idx_up].sort_values(by="auroc", ascending=False, inplace=False)
    res_dict["up"] = pd.DataFrame(df_up if head is None else df_up.iloc[0:head])

    idx_down = idx & (df["auroc"].values < 0.5)
    df_down = df.loc[idx_down].sort_values(by="auroc", ascending=True, inplace=False)
    res_dict["down"] = pd.DataFrame(df_down if head is None else df_down.iloc[0:head])



def markers(
    data: Union[MultimodalData, UnimodalData, AnnData],
    head: int = None,
    de_key: str = "de_res",
    alpha: float = 0.05,
) -> Union[Dict[str, Dict[str, pd.DataFrame]], Dict[str, Dict[str, Dict[str, pd.DataFrame]]]]:
    """ Extract DE results into a human readable structure.

    This function extracts information from ``data.varm[de_key]``, and return as
    a human readible dictionary of pandas DataFrame objects.

    Parameters
    ----------
    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``
        Data matrix with rows for cells and columns for genes.

    head: ``int``, optional, default: ``None``
        List only top ``head`` genes for each cluster. If ``None``, show any DE genes.

    de_key: ``str``, optional, default, ``de_res``
        Keyword of DE result stored in ``data.varm``.

    alpha: ``float``, optional, default: ``0.05``
        q-value threshold for getting significant DE genes. Only those with q-value of MWU test no less than ``alpha`` are significant, and thus considered as DE genes.

    Returns
    -------
    results: ``Dict[str, Dict[str, pd.DataFrame]]``
        A Python dictionary containing markers.
        If DE is performed between clusters, the structure is ``dict[cluster_id]['up' or 'down'][dataframe]``.
        If DE is performed between conditions within each cluster, the structure is ``dict[cluster_id][condition_id]['up' or 'down'][dataframe]``.
        'up' refers to up-regulated genes, which should have 'auroc' > 0.5.
        'down' refers to down-regulated genes, which should have 'auroc' < 0.5.

    Examples
    --------
    >>> marker_dict = pg.markers(data)
    """
    if de_key not in data.varm.keys():
        raise ValueError("Please run de_analysis first!")

    rec_array = data.varm[de_key]
    gene_names = data.var_names.copy()
    gene_names.name = "feature"

    de_clust = len(rec_array.dtype.names[0].split(":")) == 2

    from collections import defaultdict

    if de_clust:
        clust2cols = defaultdict(list)
        for name in rec_array.dtype.names:
            clust_id, col_name = name.split(":")
            clust2cols[clust_id].append(col_name)
        results = defaultdict(dict)
        for clust_id, col_names in clust2cols.items():
            _assemble_df(res_dict=results[clust_id],
                         rec_array=rec_array,
                         prefix=clust_id,
                         col_names=col_names,
                         gene_names=gene_names,
                         alpha=alpha,
                         head=head,
                    )
    else:
        clust2cond2cols = defaultdict(lambda: defaultdict(list))
        for name in rec_array.dtype.names:
            clust_id, cond_id, col_name = name.split(":")
            clust2cond2cols[clust_id][cond_id].append(col_name)
        results = defaultdict(lambda: defaultdict(dict))
        for clust_id, cond2cols in clust2cond2cols.items():
            for cond_id, col_names in cond2cols.items():
                _assemble_df(res_dict=results[clust_id][cond_id],
                             rec_array=rec_array,
                             prefix=":".join([clust_id, cond_id]),
                             col_names=col_names,
                             gene_names=gene_names,
                             alpha=alpha,
                             head=head,
                    )

    return results


@timer(logger=logger)
def write_results_to_excel(
    results: Union[Dict[str, Dict[str, pd.DataFrame]], Dict[str, Dict[str, Dict[str, pd.DataFrame]]]],
    output_file: str,
    ndigits: int = 3,
) -> None:
    """ Write DE analysis results into Excel workbook.

    Parameters
    ----------
    results: ``Dict[str, Dict[str, pd.DataFrame]]``, or ``Dict[str, Dict[str, Dict[str, pd.DataFrame]]]``
        DE marker dictionary generated by ``pg.markers``.

    output_file: ``str``
        File name to which the marker dictionary is written.

    ndigits: ``int``, optional, default: ``3``
        Round non p-values and q-values to ``ndigits`` after decimal point in the excel.

    Returns
    -------
    ``None``

    Marker information is written to file with name ``output_file``.

    In the generated Excel workbook,
        * If ``condition`` is ``None`` in ``pg.de_analysis``: Each tab stores DE result
          of up/down-regulated genes of cells within one cluster, and its name follows the pattern:
          **"cluster_id|up"** or **"cluster_id|dn"**.

        * If ``condition`` is not ``None`` in ``pg.de_analysis``: Each tab stores DE result
          of up/down-regulated genes of cells within one cluster under one condition level. The tab's
          name follows the pattern: **"cluster_id|cond_level|up"** or **"cluster_id|cond_level|dn"**.

        * Notice that the tab name in Excel only allows at most 31 characters. Therefore, some of the
          resulting tab names may be truncated if their names are longer than this threshold.

    Examples
    --------
    >>> pg.write_results_to_excel(marker_dict, "result.de.xlsx")
    """
    import xlsxwriter
    from natsort import natsorted

    # Need to make sure tab characters < 32 see here: https://xlsxwriter.readthedocs.io/workbook.html#:~:text=The%20worksheet%20name%20must%20be,be%20less%20than%2032%20characters.
    def compute_abbrevs(clust_ids: List[str], clust_de: bool, cond_ids: List[str], clust_abbrevs: List[str], cond_abbrevs: List[str]) -> None:
        COND_MLEN=7 # condition name should take at most 7 characters
        TOTAL_LEN=28 # excluding the last 3 |up or |dn

        clust_mlen = TOTAL_LEN
        if not clust_de:
            max_cond_len = max([len(x) for x in cond_ids])
            max_cond_len = min(max_cond_len, COND_MLEN)
            for cond_id in cond_ids:
                cond_abbrevs.append(cond_id[:max_cond_len])
            clust_mlen -= max_cond_len + 1

        for clust_id in clust_ids:
            if len(clust_id) <= clust_mlen:
                clust_abbrevs.append(clust_id)
            else:
                # Extract suffix like '-2', in which Pegasus uses to distinguish clusters with the same cell type.
                import re
                suf = ""
                match = re.search("\-\d+$", clust_id)
                if match is not None:
                    suf = match.group()
                    clust_id = clust_id[:match.start()]
                # Split by space
                fields = re.split("\s+", clust_id)
                if fields[-1].lower() in ["cell", "cells"]:
                    fields = fields[:-1]
                abbr = " ".join(fields) + suf
                if len(abbr) <= clust_mlen:
                    clust_abbrevs.append(abbr)
                else:
                    fsize = len(fields)
                    len_left = clust_mlen - len(suf) - (fsize - 1)
                    assert len_left >= fsize
                    endpos = [0] * fsize
                    while len_left > 0:
                        for i in range(fsize-1, -1, -1):
                            if endpos[i] < len(fields[i]):
                                endpos[i] += 1
                                len_left -= 1
                                if len_left == 0:
                                    break
                    clust_abbrevs.append(" ".join([fields[i][:endpos[i]] for i in range(fsize)]) + suf)

    def format_short_output_cols(
        df_orig: pd.DataFrame, ndigits: int = 3
    ) -> pd.DataFrame:
        """ Round related float columns to ndigits decimal points.
        """
        df = pd.DataFrame(df_orig, copy = True) # copy must be true, otherwise the passed df_orig will be modified.

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

    clust_ids = natsorted(results.keys())
    clust_de = isinstance(list(list(results.values())[0].values())[0], pd.DataFrame)
    cond_ids = natsorted(list(results.values())[0].keys()) if not clust_de else None
    clust_abbrevs = []
    cond_abbrevs = [] if not clust_de else None
    compute_abbrevs(clust_ids, clust_de, cond_ids, clust_abbrevs, cond_abbrevs)

    for clust_id, clust_abbr in zip(clust_ids, clust_abbrevs):
        if clust_de:
            add_worksheet(workbook, results[clust_id]["up"], f"{clust_abbr}|up")
            add_worksheet(workbook, results[clust_id]["down"], f"{clust_abbr}|dn")
        else:
            cond_dict = results[clust_id]
            for cond_id, cond_abbr in zip(cond_ids, cond_abbrevs):
                add_worksheet(workbook, cond_dict[cond_id]["up"], f"{clust_abbr}|{cond_abbr}|up")
                add_worksheet(workbook, cond_dict[cond_id]["down"], f"{clust_abbr}|{cond_abbr}|dn")

    workbook.close()
    logger.info("Excel spreadsheet is written.")


@timer(logger=logger)
def run_de_analysis(
    input_file: str,
    output_excel_file: str,
    cluster: str,
    condition: Optional[str] = None,
    de_key: str = "de_res",
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
        condition=condition,
        de_key=de_key,
        n_jobs=n_jobs,
        t=t,
        fisher=fisher,
        temp_folder=temp_folder,
        verbose=verbose,
    )

    write_output(data, input_file)
    logger.info(f"Differential expression results are written to varm/{de_key}.")

    results = markers(data, de_key=de_key, alpha=alpha)
    write_results_to_excel(results, output_excel_file, ndigits=ndigits)
