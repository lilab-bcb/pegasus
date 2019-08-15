import numpy as np
import pandas as pd
from typing import List

from scCloud.io import read_input


def search_genes(data: 'AnnData', gene_list: List[str], rec_key: str = "de_res", measure: str = "percentage") -> pd.DataFrame:
    """Extract and display gene expressions for each cluster from an `anndata` object.

    This function helps to see marker expressions in clusters via the interactive python environment.

    Parameters
    ----------

    data : `anndata` object
        An `anndata` object containing the expression matrix and differential expression results.
    gene_list : `list[str]`
        A list of gene symbols.
    rec_key : `str`
        varm keyword that stores DE results.
    measure : `str`
        Can be either `percentage` or `mean_logExpr`. `percentage` shows the percentage of cells expressed the genes and `mean_logExpr` shows the mean log expression.

    Returns
    -------
    `pandas.DataFrame`
        A data frame containing marker expressions in each cluster.

    Examples
    --------
    >>> results = misc.search_genes(data, ['CD3E', 'CD4', 'CD8'], measure = 'percentage')
    """

    columns = [x for x in data.varm[rec_key].dtype.names if x.startswith(measure + ':')]
    df = pd.DataFrame(data = data.varm[rec_key][columns], index = data.var_names)
    return df.reindex(index = gene_list)


def search_de_genes(data: 'AnnData', gene_list: List[str], rec_key: str = "de_res", test: str = "fisher", thre: float = 1.5) -> pd.DataFrame:
    """Extract and display differential expression analysis results of markers for each cluster from an `anndata` object.

    This function helps to see if markers are up or down regulated in each cluster via the interactive python environment. `++` indicates up-regulated and fold change >= threshold, `+` indicates up-regulated but fold change < threshold, `--` indicates down-regulated and fold change <= 1 / threshold, `-` indicates down-regulated but fold change > 1 / threshold, '?' indicates not differentially expressed.

    Parameters
    ----------

    data : `anndata` object
        An `anndata` object containing the expression matrix and differential expression results.
    gene_list : `list[str]`
        A list of gene symbols.
    rec_key : `str`
        varm keyword that stores DE results.
    test : `str`, optional (default: `fisher`)
        Differential expression test to look at, could be either `t`, `fisher` or `mwu`.
    thre : `float`, optional (default: `1.5`)
        Fold change threshold to determine if the marker is a strong DE (`++` or `--`) or weak DE (`+` or `-`).

    Returns
    -------
    `pandas.DataFrame`
        A data frame containing marker differential expression results for each cluster.

    Examples
    --------
    >>> results = misc.search_de_genes(data, ['CD3E', 'CD4', 'CD8'], test = 'fisher', thre = 2.0)
    """

    columns = [x for x in data.varm[rec_key].dtype.names if x.startswith(test + '_qval:')]
    df_de = pd.DataFrame(data.varm[rec_key][columns], index = data.var_names)
    df_de = df_de.reindex(index = gene_list)

    columns = [x for x in data.varm[rec_key].dtype.names if (x.startswith("percentage_fold_change:") if test == "fisher" else x.startswith("log_fold_change:"))]    
    df_fc = pd.DataFrame(data.varm[rec_key][columns], index = data.var_names)
    df_fc = df_fc.reindex(index = gene_list)
    if test != "fisher":
        df_fc = np.exp(df_fc)

    results = np.zeros((len(gene_list), len(columns)), dtype=np.dtype("U4"))
    results[:] = "?"
    results[np.isnan(df_de)] = "NaN"
    results[(df_de <= 0.05).values & (df_fc > 1.0).values] = "+"
    results[(df_de <= 0.05).values & (df_fc >= thre).values] = "++"
    results[(df_de <= 0.05).values & (df_fc < 1.0).values] = "-"
    results[(df_de <= 0.05).values & (df_fc <= 1.0 / thre).values] = "--"

    clusts = [x.rpartition(':')[2] for x in columns]
    df = pd.DataFrame(data=results, index=gene_list, columns=clusts)
    return df


def show_attributes(
    input_file: str, show_attributes: bool, show_gene_attributes: bool, show_values_for_attributes: str
) -> None:
    """ Show data attributes. For command line use.
    """

    data = read_input(input_file, h5ad_mode="r")
    if show_attributes:
        print(
            "Available sample attributes in input dataset: {0}".format(
                ", ".join(data.obs.columns.values)
            )
        )
    if show_gene_attributes:
        print(
            "Available gene attributes in input dataset: {0}".format(
                ", ".join(data.var.columns.values)
            )
        )
    if not show_values_for_attributes is None:
        for attr in show_values_for_attributes.split(","):
            print(
                "Available values for attribute {0}: {1}.".format(
                    attr, ", ".join(np.unique(data.obs[attr]))
                )
            )


def perform_oneway_anova(data, glist, restriction_vec, group_str, fdr_alpha=0.05):
    from scipy.stats import f_oneway
    from statsmodels.stats.multitest import fdrcorrection as fdr

    selected = np.ones(data.shape[0], dtype=bool)
    for rest_str in restriction_vec:
        attr, value_str = rest_str.split(":")
        values = value_str.split(",")
        selected = selected & np.isin(data.obs[attr], values)
    gene_list = np.array(glist)
    gene_list = gene_list[np.isin(gene_list, data.var_names)]
    newdat = data[selected, :][:, gene_list].copy()
    newdat.X = newdat.X.toarray()
    group_attr, tmp_str = group_str.split(":")
    groups_str = tmp_str.split(";")
    ngr = len(groups_str)
    group_names = []
    group_idx = np.zeros((ngr, newdat.shape[0]), dtype=bool)
    for i, gstr in enumerate(groups_str):
        name, values = gstr.split("~")
        group_names.extend([name + "_mean", name + "_percent"])
        group_idx[i] = np.isin(newdat.obs[group_attr], values.split(","))
    np.warnings.filterwarnings("ignore")
    stats = np.zeros((len(gene_list), 3 + ngr * 2))
    for i in range(len(gene_list)):
        arr_list = []
        for j in range(group_idx.shape[0]):
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
    cols.extend(group_names)
    raw_results = pd.DataFrame(stats, columns=cols, index=gene_list)
    results = raw_results[raw_results["qval"] <= fdr_alpha]
    results = results.sort_values("qval")
    return results, raw_results
