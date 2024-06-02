import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd

from pegasus.tools import predefined_pathways, load_signatures_from_file, eff_n_jobs
from pegasusio import MultimodalData, UnimodalData, timer
from typing import Union, Optional, Tuple


@timer(logger=logger)
def gsea(
    data: Union[MultimodalData, UnimodalData],
    de_key: str,
    pathways: str,
    method: str = "blitzgsea",
    gsea_key: str = "gsea_out",
    min_size: Optional[int] = 5,
    max_size: Optional[int] = 4000,
    n_jobs: int = 1,
    seed: int = 0,
    **kwargs,
) -> None:
    """Perform Gene Set Enrichment Analysis (GSEA).

    Parameters
    ----------
    data: Union[``MultimodalData``, ``UnimodalData``]
        Single-cell or pseudo-bulk data.

    de_key: ``str``
        Key in pre-computed DE results representing gene ranks.
        The string format is ``res_key:attr``, where ``res_key`` is the key of DE results in ``data.varm``, and ``attr`` is the column name in ``data.varm['res_key']`` used as the gene signatures for GSEA analysis.

    pathways: ``str``
        Either a string or a path to the gene set file in GMT format. If string, choosing from "hallmark" and "canonical_pathways" (MSigDB H and C2/CP).

    method: ``str``, optional, default: ``blitzgsea``
        Specify which package to use as the backend for GSEA.
        By default, use ``blitzgsea`` which is a pure Python fast implementation of the pre-rank GSEA algorithm.
        Alternatively, if specify ``fgsea``, then use R package ``fgsea``, which requires ``rpy2`` and R installation.

    gsea_key: ``str``, optional, default: ``"gsea_out"``
        Key to use to store GSEA results as a data frame.

    min_size: ``int``, optional, default: ``5``
        Minimal size of a gene set to consider.

    max_size: ``int``, optional, default: ``4000``
        Maximal size of a gene set to consider.

    n_jobs: ``int``, optional, default: ``1``
        Numbr of processes for parallel computation.

    seed: ``int``, optional, default: ``0``
        Random seed to make sure GSEA results are reproducible.

    Returns
    -------
    ``None``

    Update ``data.uns``:
        ``data.uns[gsea_key]``: GSEA outputs sorted by adjusted p-values.

    Examples
    --------
    >>> pg.gsea(data, de_key='deseq2:stat', pathways='hallmark')
    >>> pg.gsea(data, de_key='de_res:1:mwu_U', pathways="custom.gmt", method="fgsea")
    """
    if method == "blitzgsea":
        # center default is False, unless explicitly specified by users
        permutations = 2000
        center = False
        if 'permutations' in kwargs:
            permutations = kwargs['permutations']
            del kwargs['permutations']
        if 'center' in kwargs:
            center = kwargs['center']
            del kwargs['center']

        _run_blitzgsea(
            data=data,
            de_key=de_key,
            pathways=pathways,
            gsea_key=gsea_key,
            min_size=min_size,
            max_size=max_size,
            n_jobs=n_jobs,
            seed=seed,
            permutations=permutations,
            center=center,
            **kwargs,
        )
    else:
        _run_fgsea(
            data=data,
            de_key=de_key,
            pathways=pathways,
            gsea_key=gsea_key,
            min_size=min_size,
            max_size=max_size,
            n_jobs=n_jobs,
            seed=seed,
        )


def _split_de_key(de_key: str) -> Tuple[str, str, Optional[str]]:
    assert ":" in de_key, f"Key '{de_key}' does not satisfy format 'df_key:attr'!"
    keys = de_key.split(":")
    df_key = keys[0]
    rank_key = ":".join(keys[1:])
    cluster_key = None if len(keys) <= 2 else keys[1]

    return (df_key, rank_key, cluster_key)


def _run_blitzgsea(
    data: Union[MultimodalData, UnimodalData],
    de_key: str,
    pathways: str,
    gsea_key: str,
    min_size: int,
    max_size: int,
    n_jobs: int,
    seed: int,
    permutations: int,
    center: bool,
    **kwargs,
) -> None:
    try:
        import blitzgsea as blitz
    except ModuleNotFoundError as e:
        import sys
        logger.error(f"{e}\nNeed blitzgsea! Try 'pip install blitzgsea'.")
        sys.exit(-1)

    df_key, rank_key, cluster_key = _split_de_key(de_key)
    assert df_key in data.varm, f"Key '{df_key}' not in data.varm! Wrong key name or need to run DE analysis first!"
    assert rank_key in data.varm[df_key].dtype.names, f"Key '{rank_key}' not in DE result! Wrong key name specified!"

    # Decide q-value column name
    qstr = "padj"
    if cluster_key:
        # Not pseudobulk DE
        test_str = rank_key.split(":")[-1].split("_")[0]
        qstr = f"{test_str}_qval"
        assert qstr in data.varm[df_key].dtype.names, f"Q-value key '{qstr}' not in DE result! Either q-values are missing, or need to rename the column to this name!"

    rank_df = pd.DataFrame(data.varm[df_key], index=data.var_names)
    rank_df = rank_df.loc[~np.isnan(rank_df[qstr])].copy()    # Ignore genes of NaN q-values, which is the case for Wald test in DESeq2 model.
    rank_df = rank_df[[rank_key]].reset_index()
    rank_df.columns = [0, 1]
    rank_df = rank_df.sort_values(by=1, ascending=False)

    library = load_signatures_from_file(predefined_pathways.get(pathways, pathways))

    n_jobs = eff_n_jobs(n_jobs)

    result_df = blitz.gsea(
        signature=rank_df,
        library=library,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        processes=n_jobs,
        permutations=permutations,
        center=center,
        **kwargs,
    )
    result_df.sort_values('fdr', inplace=True)
    data.uns[gsea_key] = result_df


def _run_fgsea(
    data: Union[MultimodalData, UnimodalData],
    de_key: str,
    pathways: str,
    gsea_key: str,
    min_size: str,
    max_size: str,
    n_jobs: int,
    seed: int,
):
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
        nproc = 0    # Avoid setting BPPARM when only using 1 thread

    ro.r(f"set.seed({seed})")
    pwdict = load_signatures_from_file(predefined_pathways.get(pathways, pathways))
    pathways_r = ro.ListVector(pwdict)

    df_key, rank_key, cluster_key = _split_de_key(de_key)
    assert df_key in data.varm, f"Key '{df_key}' not in data.varm! Wrong key name or run DE analysis beforehand!"
    assert rank_key in data.varm[df_key].dtype.names, f"Key '{rank_key}' not in DE result! Wrong key name specified!"

    # Decide q-value column name
    qstr = "padj"
    if cluster_key:
        # Not pseudobulk DE
        test_str = rank_key.split(":")[-1].split("_")[0]
        qstr = f"{test_str}_qval"
        assert qstr in data.varm[df_key].dtype.names, f"Q-value key '{qstr}' not in DE result! Either q-values are missing, or need to rename the column to this name!"

    qvals = data.varm[df_key][qstr]
    idx_select = np.where(~np.isnan(qvals))[0]    # Ignore genes with NaN q-values, which is the case for Wald test in DESeq2 model.
    rank_vec = ro.FloatVector(data.varm[df_key][rank_key][idx_select])
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
    # Make result data frame column names consistent with BlitzGSEA
    res_df.rename(columns={
        "pathway": "Term",
        "padj": "fdr",
        "ES": "es",
        "NES": "nes",
        "size": "geneset_size",
        "leadingEdge": "leading_edge",
    }, inplace=True)
    res_df.sort_values("fdr", inplace=True)
    res_df = res_df.set_index("Term")
    data.uns[gsea_key] = res_df


@timer(logger=logger)
def write_gsea_results_to_excel(
    data: Union[MultimodalData, UnimodalData],
    output_file: str,
    gsea_key: Optional[str] = "gsea_out",
    ndigits: Optional[int] = 3,
) -> None:
    """Write Gene Set Enrichment Analysis (GSEA) results generated by gsea function into Excel workbook.

    Parameters
    ----------
    data: Union[``MultomodalData``, ``UnimodalData``]
        Single-cell or pseudo-bulk data.

    output_file: ``str``
        File name for the output.

    gsea_key: ``str``, optinoal, default: ``"gsea_out"``
        Key name of GSEA results stored in ``data.uns`` field.

    ndigits: ``int``, optional, default: ``3``
        Round non p-values and q-values to ``ndigits`` after decimal point in the excel.

    Returns
    -------
    ``None``

    Pathway information is written to file with name ``output_file``.

    In the generated Excel workbook, there are two tabs:
        * ``UP``: Pathways with positive NES scores (``nes`` column), which are sorted by q-values (``fdr`` column).
        * ``DOWN``: Pathways with negative NES scores (``nes`` column), which are sorted by q-values (``fdr`` column).

    Examples
    --------
    >>> pg.write_gsea_results_to_excel(data, "gsea.xlsx")
    """
    assert gsea_key in data.uns.keys(), f"Key '{gsea_key}' does not exist in data.uns!"

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

        df_gsea = pd.DataFrame(data.uns[gsea_key])
        add_worksheet(workbook, df_gsea.loc[df_gsea['nes']>0].set_index("Term"), "UP")
        add_worksheet(workbook, df_gsea.loc[df_gsea['nes']<0].set_index("Term"), "DOWN")
        logger.info("Excel spreadsheet is written.")
    finally:
        workbook.close()
