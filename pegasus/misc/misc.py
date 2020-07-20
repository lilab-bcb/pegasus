import numpy as np
import pandas as pd
from typing import List
from anndata import AnnData

import logging
logger = logging.getLogger("pegasus")



def search_genes(
    data: AnnData,
    gene_list: List[str],
    rec_key: str = "de_res",
    measure: str = "percentage",
) -> pd.DataFrame:
    """Extract and display gene expressions for each cluster from an `anndata` object.

    This function helps to see marker expressions in clusters via the interactive python environment.

    Parameters
    ----------

    data: ``anndata.AnnData``
        Annotated data matrix containing the expression matrix and differential expression results.

    gene_list: ``List[str]``
        A list of gene symbols.

    rec_key: ``str``, optional, default: ``"de_res"``
        Keyword of DE analysis result stored in ``data.varm``.

    measure : ``str``, optional, default: ``"percentage"``
        Can be either ``"percentage"`` or ``"mean_logExpr"``:
            * ``percentage`` shows the percentage of cells expressed the genes;
            * ``mean_logExpr`` shows the mean log expression.

    Returns
    -------
    ``pandas.DataFrame``
        A data frame containing marker expressions in each cluster.

    Examples
    --------
    >>> results = pg.search_genes(adata, ['CD3E', 'CD4', 'CD8'])
    """

    columns = [x for x in data.varm[rec_key].dtype.names if x.startswith(measure + ":")]
    df = pd.DataFrame(data=data.varm[rec_key][columns], index=data.var_names)
    return df.reindex(index=gene_list)


def search_de_genes(
    data: AnnData,
    gene_list: List[str],
    rec_key: str = "de_res",
    de_test: str = "fisher",
    de_alpha: float = 0.05,
    thre: float = 1.5,
) -> pd.DataFrame:
    """Extract and display differential expression analysis results of markers for each cluster.

    This function helps to see if markers are up or down regulated in each cluster via the interactive python environment:
        * ``++`` indicates up-regulated and fold change >= threshold;
        * ``+`` indicates up-regulated but fold change < threshold;
        * ``--`` indicates down-regulated and fold change <= 1 / threshold;
        * ``-`` indicates down-regulated but fold change > 1 / threshold;
        * ``?`` indicates not differentially expressed.

    Parameters
    ----------
    data: ``anndata.Anndata``
        Annotated data matrix containing the expression matrix and differential expression results.

    gene_list: ``List[str]``
        A list of gene symbols.

    rec_key: ``str``, optional, default: ``"de_res"``
        Keyword of DE analysis result stored in ``data.varm``.

    de_test : ``str``, optional, default: ``"fisher"``
        Differential expression test to look at, could be either ``t``, ``fisher`` or ``mwu``.

    de_alpha : ``float``, optional, default: ``0.05``
        False discovery rate.

    thre : ``float``, optional, default: ``1.5``
        Fold change threshold to determine if the marker is a strong DE (``++`` or ``--``) or weak DE (``+`` or ``-``).

    Returns
    -------
    ``pandas.DataFrame``
        A data frame containing marker differential expression results for each cluster.

    Examples
    --------
    >>> df = pegasus.misc.search_de_genes(adata, ['CD3E', 'CD4', 'CD8'], thre = 2.0)
    """

    columns = [
        x for x in data.varm[rec_key].dtype.names if x.startswith(de_test + "_qval:")
    ]
    df_de = pd.DataFrame(data.varm[rec_key][columns], index=data.var_names)
    df_de = df_de.reindex(index=gene_list)

    columns = [
        x
        for x in data.varm[rec_key].dtype.names
        if (
            x.startswith("percentage_fold_change:")
            if de_test == "fisher"
            else x.startswith("log_fold_change:")
        )
    ]
    df_fc = pd.DataFrame(data.varm[rec_key][columns], index=data.var_names)
    df_fc = df_fc.reindex(index=gene_list)
    if de_test != "fisher":
        df_fc = np.exp(df_fc)

    results = np.zeros((len(gene_list), len(columns)), dtype=np.dtype("U4"))
    results[:] = "?"
    results[np.isnan(df_de)] = "NaN"
    results[(df_de <= de_alpha).values & (df_fc > 1.0).values] = "+"
    results[(df_de <= de_alpha).values & (df_fc >= thre).values] = "++"
    results[(df_de <= de_alpha).values & (df_fc < 1.0).values] = "-"
    results[(df_de <= de_alpha).values & (df_fc <= 1.0 / thre).values] = "--"

    clusts = [x.rpartition(":")[2] for x in columns]
    df = pd.DataFrame(data=results, index=gene_list, columns=clusts)
    return df


def show_attributes(
    input_file: str,
    show_attributes: bool,
    show_gene_attributes: bool,
    show_values_for_attributes: str,
) -> None:
    """ Show data attributes. For command line use.
    """

    # data = read_input(input_file, mode="r")
    # if show_attributes:
    #     print(
    #         "Available sample attributes in input dataset: {0}".format(
    #             ", ".join(data.obs.columns.values)
    #         )
    #     )
    # if show_gene_attributes:
    #     print(
    #         "Available gene attributes in input dataset: {0}".format(
    #             ", ".join(data.var.columns.values)
    #         )
    #     )
    # if not show_values_for_attributes is None:
    #     for attr in show_values_for_attributes.split(","):
    #         print(
    #             "Available values for attribute {0}: {1}.".format(
    #                 attr, ", ".join(np.unique(data.obs[attr]))
    #             )
    #         )


def perform_oneway_anova(
    data: AnnData,
    glist: List[str],
    restriction_vec: List[str],
    group_str: str,
    fdr_alpha: float = 0.05,
    res_key: str = None,
) -> pd.DataFrame:
    """Perform one way ANOVA on a subset of cells (restricted by restriction_vec) grouped by group_str and control FDR at fdr_alpha.
    Parameters
    ----------

    data : `anndata` object
        An `anndata` object containing the expression matrix.
    glist : `list[str]`
        A list of gene symbols.
    restriction_vec : `list[str]`
        A vector of restrictions for selecting cells. Each restriction takes the format of attr:value,value,value
    group_str : `str`
        How to group selected cells for ANOVA analysis. If group_str is for pseudotime, it has two formats. 1) 'pseudotime:time:n', which divides cells by equal pseudotime invertal; 2) 'pseudotime:size:n' divides cells by equal number of cells.
    fdr_alpha : `float`, optional (default: 0.05)
        False discovery rate.
    res_key : `str`, optional (default: None)
        Store results into data using res_key, the grouping information is stored in obs and the results is stored in uns.

    Returns
    -------
    `pandas.DataFrame`
        Results for genes that pass FDR control.

    Examples
    --------
    >>> results = misc.perform_oneway_anova(data, ['CD3E', 'CD4', 'CD8'], [], 'pseudotime:size:10')
    """

    from scipy.stats import f_oneway
    from statsmodels.stats.multitest import fdrcorrection as fdr

    selected = np.ones(data.shape[0], dtype=bool)
    for rest_str in restriction_vec:
        attr, value_str = rest_str.split(":")
        values = value_str.split(",")
        selected = selected & np.isin(data.obs[attr], values)

    gene_list = np.array(glist)
    gene_list = gene_list[np.isin(gene_list, data.var_names)]
    ngene = gene_list.size

    newdat = data[selected, :][:, gene_list].copy()
    newdat.X = newdat.X.toarray()

    group_values = group_str.split(":")
    group_names = []
    col_names = []

    ngr = 0
    group_idx = None

    if group_values[0] == "pseudotime":
        assert len(group_values) == 3
        div_by = group_values[1]
        ngr = int(group_values[2])

        group_idx = np.zeros((ngr, newdat.shape[0]), dtype=bool)
        pseudotimes = newdat.obs["pseudotime"].values

        min_t = pseudotimes.min()
        max_t = pseudotimes.max()

        if div_by == "time":
            interval = (max_t - min_t) / ngr
            left = min_t - 1e-5
            for i in range(ngr):
                right = min_t + interval * (i + 1)
                name = "({:.2f}, {:.2f}]".format(left if left >= 0 else 0.0, right)
                group_names.append(name)
                group_idx[i] = (pseudotimes > left) & (pseudotimes <= right)
                left = right
        else:
            assert div_by == "size"
            ords = np.argsort(pseudotimes)
            quotient = ords.size // ngr
            residule = ords.size % ngr

            fr = 0
            for i in range(ngr):
                to = fr + quotient + (i < residule)
                name = "[{:.2f}, {:.2f}]".format(
                    pseudotimes[ords[fr]], pseudotimes[ords[to - 1]]
                )
                group_names.append(name)
                group_idx[i][ords[fr:to]] = True
                fr = to

    else:
        assert len(group_values) == 2
        group_attr = group_values[0]
        tmp_str = group_values[1]
        groups_str = tmp_str.split(";")

        ngr = len(groups_str)
        group_idx = np.zeros((ngr, newdat.shape[0]), dtype=bool)

        for i, gstr in enumerate(groups_str):
            name, values = gstr.split("~")
            group_names.append(name)
            group_idx[i] = np.isin(newdat.obs[group_attr], values.split(","))

    for i in range(ngr):
        print("Group {} has {} cells.".format(group_names[i], group_idx[i].sum()))

    np.warnings.filterwarnings("ignore")
    stats = np.zeros((ngene, 3 + ngr * 2))
    for i in range(ngene):
        arr_list = []
        for j in range(ngr):
            arr = newdat.X[group_idx[j], i]
            stats[i, 3 + j * 2] = arr.mean()
            stats[i, 3 + j * 2 + 1] = (arr > 0).sum() * 100.0 / arr.size
            arr_list.append(arr)
        stats[i, 0], stats[i, 1] = f_oneway(*arr_list)
        if np.isnan(stats[i, 0]):
            stats[i, 0] = 0.0
            stats[i, 1] = 1.0
    passed, stats[:, 2] = fdr(stats[:, 1])

    cols = ["fstat", "pval", "qval"]
    for i in range(ngr):
        cols.extend([group_names[i] + "_mean", group_names[i] + "_percent"])
    raw_results = pd.DataFrame(stats, columns=cols, index=gene_list)

    results = raw_results[raw_results["qval"] <= fdr_alpha]
    results = results.sort_values("qval")

    if res_key is not None:
        data.uns[res_key] = raw_results
        data.obs[res_key] = "background"
        for i in range(ngr):
            idx = np.zeros(data.shape[0], dtype=bool)
            idx[selected] = group_idx[i]
            data.obs.loc[idx, res_key] = group_names[i]

    return results
