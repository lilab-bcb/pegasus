import logging

logger = logging.getLogger(__name__)

import pandas as pd
import numpy as np

from pegasus.tools import predefined_pathways, load_signatures_from_file, eff_n_jobs
from pegasusio import MultimodalData, UnimodalData, timer
from typing import Union, Optional


@timer(logger=logger)
def gsea(
    data: Union[MultimodalData, UnimodalData],
    rank_key: str,
    pathways: str,
    de_key: str = "de_res",
    method: str = "gseapy",
    gsea_key: str = "gsea_out",
    min_size: int = 15,
    max_size: int = 500,
    n_jobs: int = 4,
    seed: int = 0,
    verbose: bool = True,
    **kwargs,
) -> None:
    """Perform Gene Set Enrichment Analysis (GSEA).

    Parameters
    ----------
    data: ``MultimodalData`` or ``UnimodalData``
        Single-cell or pseudo-bulk data.

    rank_key: ``str``
        Key in pre-computed DE results representing gene ranks.
        The string format is ``de_key:attr``, where ``de_key`` is the key of DE results in ``data.varm``, and ``attr`` is the column name in ``data.varm['de_key']`` used as the gene signatures for GSEA analysis.

    pathways: ``str``
        Either a keyword or a path to the gene set file in GMT format. If keyword, choosing from "hallmark" and "canonical_pathways" (MSigDB H and C2/CP).

    de_key: ``str``, optional, default: ``de_res``
        Key name of DE analysis results stored.`data.varm[de_key]` should contain a record array of DE results.

    method: ``str``, optional, default: ``gseapy``
        Specify which package to use as the backend for GSEA.
        By default ``gseapy``, use GSEAPY's prerank method. Notice that ``permutation_num=1000`` is set by default. If you want to change these parameters, please reset in ``kwargs``.
        Alternatively, if specify ``fgsea``, then use R package ``fgsea``, which requires ``rpy2`` and R installation.

    gsea_key: ``str``, optional, default: ``"gsea_out"``
        Key to use to store GSEA results as a data frame.

    min_size: ``int``, optional, default: ``15``
        Minimum allowed number of genes from gene set also the data set.

    max_size: ``int``, optional, default: ``500``
        Maximum allowed number of genes from gene set also the data set.

    n_jobs: ``int``, optional, default: ``4``
        Numbr of threads used for parallel computation.

    seed: ``int``, optional, default: ``0``
        Random seed to make sure GSEA results are reproducible.

    verbose: ``bool``, optional, default: ``True``
        If printing out progress of the job. Only works when ``method="gseapy"``.

    kwargs
        If ``method="gseapy"``, pass other keyword arguments to ``gseapy.prerank`` function. Details about GSEAPY prerank function's optional parameters are `here <https://gseapy.readthedocs.io/en/latest/run.html#gseapy.prerank>`_.

    Returns
    -------
    ``None``

    Update ``data.uns``:
        ``data.uns[gsea_key]``: GSEA outputs sorted by adjusted p-values.

    Examples
    --------
    >>> pg.gsea(data, "deseq2:stat", "canonical_pathways")
    >>> pg.gsea(data, "de_res:1:mwu_U", "canonical_pathways", method="fgsea")
    """
    if method == "gseapy":
        _run_gseapy(
            data=data,
            rank_key=rank_key,
            pathways=pathways,
            de_key=de_key,
            gsea_key=gsea_key,
            min_size=min_size,
            max_size=max_size,
            n_jobs=n_jobs,
            seed=seed,
            verbose=verbose,
            **kwargs,
        )
    else:
        _run_fgsea(
            data=data,
            rank_key=rank_key,
            pathways=pathways,
            de_key=de_key,
            gsea_key=gsea_key,
            min_size=min_size,
            max_size=max_size,
            n_jobs=n_jobs,
            seed=seed,
        )


def _decide_qval_str(rank_key):
    qstr = "padj"
    if ":" in rank_key:
        # pegasus.de_analysis result
        prefix = rank_key.split("_")[0]
        qstr = f"{prefix}_qval"
    return qstr


def _validate_keys(data, de_key, rank_key, padj_key):
    if de_key not in data.varm:
        import sys
        logger.error(f"Key '{de_key}' not in data.varm! Wrong key name or need to run DE analysis first!")
        sys.exit(-1)
    if rank_key not in data.varm[de_key].dtype.names:
        import sys
        logger.error(f"Key '{rank_key}' not in DE result! Wrong key name specified!")
        sys.exit(-1)
    if padj_key not in data.varm[de_key].dtype.names:
        import sys
        logger.error(f"Q-value key '{padj_key}' not in DE result! Either q-values are missing, or need to rename the column to this name!")
        sys.exit(-1)


def _run_gseapy(
    data: Union[MultimodalData, UnimodalData],
    rank_key: str,
    pathways: str,
    de_key: str,
    gsea_key: str,
    min_size: int,
    max_size: int,
    n_jobs: int,
    seed: int,
    verbose: bool,
    **kwargs,
) -> None:
    try:
        import gseapy as gp
    except ModuleNotFoundError as e:
        import sys
        logger.error(f"{e}\nNeed gseapy! Try 'pip install gseapy'.")
        sys.exit(-1)

    qstr = _decide_qval_str(rank_key)
    _validate_keys(data, de_key, rank_key, qstr)

    qvals = data.varm[de_key][qstr]
    idx_select = np.where(~np.isnan(qvals))[0]    # Ignore genes with NaN q-values, which is the case for independent filtering in DESeq2 model.
    rank_df = pd.DataFrame({'gene': data.var_names[idx_select], 'rank': data.varm[de_key][rank_key][idx_select]})

    gene_sets = load_signatures_from_file(predefined_pathways.get(pathways, pathways))
    n_jobs = eff_n_jobs(n_jobs)
    res = gp.prerank(
        rnk=rank_df,
        gene_sets=gene_sets,
        min_size=min_size,
        max_size=max_size,
        threads=n_jobs,
        seed=seed,
        verbose=verbose,
        **kwargs,
    )

    res_df = res.res2d
    res_df.rename(columns={"FDR q-val": "padj", "Term": "pathway"}, inplace=True)
    res_df["NES"] = res_df["NES"].astype(np.float64)
    res_df["padj"] = res_df["padj"].astype(np.float64)
    res_df.sort_values(["padj", "pathway"], ascending=[True, True], inplace=True)
    data.uns[gsea_key] = res_df


def _run_fgsea(
    data: MultimodalData,
    rank_key: str,
    pathways: str,
    de_key: str,
    gsea_key: str,
    min_size: int,
    max_size: int,
    n_jobs: int,
    seed: int,
) -> None:
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter
    except ModuleNotFoundError as e:
        import sys

        logger.error(f"{e}\nNeed rpy2! Try 'pip install rpy2'.")
        sys.exit(-1)

    try:
        fgsea = importr("fgsea")
    except ModuleNotFoundError:
        import sys

        text = """Please install fgsea in order to run this function.\n
                To install this package, start R and enter:\n
                if (!require("BiocManager", quietly = TRUE))
                    install.packages("BiocManager")
                BiocManager::install("fgsea")"""

        logger.error(text)
        sys.exit(-1)

    nproc = eff_n_jobs(n_jobs)
    if nproc == 1:
        nproc = 0    # Avoid setting BPPARAM when only using 1 thread

    ro.r(f"set.seed({seed})")
    pwdict = load_signatures_from_file(predefined_pathways.get(pathways, pathways))
    pathways_r = ro.ListVector(pwdict)

    qstr = _decide_qval_str(rank_key)
    _validate_keys(data, de_key, rank_key, qstr)

    qvals = data.varm[de_key][qstr]
    idx_select = np.where(~np.isnan(qvals))[0]    # Ignore genes with NaN q-values, which is the case for independent filtering in DESeq2 model.
    rank_vec = ro.FloatVector(data.varm[de_key][rank_key][idx_select])
    rank_vec.names = ro.StrVector(data.var_names[idx_select])
    res = fgsea.fgsea(pathways_r, rank_vec, minSize=min_size, maxSize=max_size, nproc=nproc)
    unlist = ro.r(
        """
        function(df) {
            df$leadingEdge <- sapply(df$leadingEdge, function(x) {paste(unlist(x), collapse=',')})
            return(df)
        }
        """
    )
    with localconverter(ro.default_converter + pandas2ri.converter):
        res_df = ro.conversion.rpy2py(unlist(res))
    res_df.sort_values(["padj", "pathway"], ascending=[True, True], inplace=True)
    data.uns[gsea_key] = res_df


@timer(logger=logger)
def write_gsea_results_to_excel(
    data: Union[MultimodalData, UnimodalData],
    output_file: str,
    gsea_key: str = "gsea_out",
    ndigits: int = 3,
) -> None:
    """Write Gene Set Enrichment Analysis (GSEA) results into Excel workbook.

    Parameters
    ----------
    data: ``MultomodalData`` or ``UnimodalData``
        Single-cell or pseudo-bulk data.

    output_file: ``str``
        File name for the output.

    gsea_key: ``str``, optinoal, default: ``gsea_out``
        Key name of GSEA results stored in ``data.uns`` field.

    ndigits: ``int``, optional, default: ``3``
        Round non p-values and q-values to ``ndigits`` after decimal point in the excel.
    """
    assert fgsea_key in data.uns.keys(), f"Key '{fgsea_key}' does not exist in data.uns!"

    import xlsxwriter

    def format_short_output_cols(
        df_orig: pd.DataFrame, ndigits: int
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
        df = format_short_output_cols(df_orig, ndigits)
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

    try:
        workbook = xlsxwriter.Workbook(output_file, {"nan_inf_to_errors": True})
        workbook.formats[0].set_font_size(9)

        df_fgsea = pd.DataFrame(data.uns[fgsea_key])
        add_worksheet(workbook, df_fgsea.loc[df_fgsea['NES']>0].set_index("pathway"), "UP")
        add_worksheet(workbook, df_fgsea.loc[df_fgsea['NES']<0].set_index("pathway"), "DOWN")
        logger.info("Excel spreadsheet is written.")
    finally:
        workbook.close()