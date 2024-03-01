import logging

logger = logging.getLogger(__name__)

import pandas as pd

from pegasus.tools import predefined_pathways, load_signatures_from_file
from pegasusio import MultimodalData, UnimodalData, timer
from typing import Union, Optional


@timer(logger=logger)
def fgsea(
    data: Union[MultimodalData, UnimodalData],
    log2fc_key: str,
    pathways: str,
    de_key: Optional[str] = "de_res",
    minSize: Optional[int] = 15,
    maxSize: Optional[int] = 500,
    nproc: Optional[int] = 0,
    seed: Optional[int] = 0,
    fgsea_key: Optional[str] = "fgsea_out",
) -> None:
    """Perform Gene Set Enrichment Analysis using fGSEA. This function calls R package fGSEA, requiring fGSEA in R installed.

    Parameters
    ----------
    data: Union[``MultimodalData``, ``UnimodalData``]
        Single-cell or pseudo-bulk data.

    log2fc_key: ``str``
        Key in pre-computed DE results representing log2 fold change.

    pathways: ``str``
        Either a string or a path to the gene set file in GMT format. If string, choosing from "hallmark" and "canonical_pathways" (MSigDB H and C2/CP).

    de_key: ``str``, optional, default: ``"de_res"``
        Key name of DE analysis results stored. data.varm[de_key] should contain a record array of DE results.

    minSize: ``int``, optional, default: ``15``
        Minimal size of a gene set to consider.

    maxSize: ``int``, optional, default: ``500``
        Maximal size of a gene set to consider.

    nproc: ``int``, optional, default: ``0``
        Numbr of processes for parallel computation. If nproc > 0, set BPPARAM.

    seed: ``int``, optional, default: ``0``
        Random seed to make sure fGSEA results are reproducible.

    fgsea_key: ``str``, optional, default: ``"fgsea_out"``
        Key to use to store fGSEA results as a data frame.

    Returns
    -------
    ``None``

    Update ``data.uns``:
        ``data.uns[fgsea_key]``: fGSEA outputs sorted by padj.

    Examples
    --------
    >>> pg.fgsea(data, '3:log2FC', hallmark', fgsea_key='fgsea_res')
    """
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

    ro.r(f"set.seed({seed})")
    pwdict = load_signatures_from_file(predefined_pathways.get(pathways, pathways))
    pathways_r = ro.ListVector(pwdict)
    log2fc = ro.FloatVector(data.varm[de_key][log2fc_key])
    log2fc.names = ro.StrVector(data.var_names)
    res = fgsea.fgsea(pathways_r, log2fc, minSize=minSize, maxSize=maxSize, nproc=nproc)
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
    res_df.sort_values("padj", inplace=True)
    data.uns[fgsea_key] = res_df


@timer(logger=logger)
def write_fgsea_results_to_excel(
    data: Union[MultimodalData, UnimodalData],
    output_file: str,
    fgsea_key: Optional[str] = "fgsea_out",
    ndigits: Optional[int] = 3,
) -> None:
    """
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