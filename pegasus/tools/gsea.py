import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd

from pegasus.tools import predefined_pathways, load_signatures_from_file, eff_n_jobs
from pegasusio import MultimodalData, UnimodalData, timer
from typing import Union, Optional


@timer(logger=logger)
def gsea(
    data: Union[MultimodalData, UnimodalData],
    pathways: str,
    de_key: str = "de_res:stat",
    method: str = "blitzgsea",
    gsea_key: str = "gsea_out",
    min_size: int = 5,
    max_size: int = 4000,
    n_jobs: int = 1,
    seed: int = 0,
    **kwargs,
) -> None:
    """Perform Gene Set Enrichment Analysis (GSEA).

    Parameters
    ----------
    data: Union[``MultimodalData``, ``UnimodalData``]
        Single-cell or pseudo-bulk data.

    rank_key: ``str``, optional, default: ``stat``
        Key in pre-computed DE results representing gene ranks.
        By default, use the test statistics in the pre-computed DE results.

    pathways: ``str``
        Either a string or a path to the gene set file in GMT format. If string, choosing from "hallmark" and "canonical_pathways" (MSigDB H and C2/CP).

    backend: ``str``, optional, default: ``blitzgsea``
        Specify which package to use as the backend for GSEA.
        By default, use ``blitzgsea`` which is a pure Python fast implementation of the pre-rank GSEA algorithm.
        Alternatively, if specify ``fgsea``, then use R package ``fgsea``, which requires ``rpy2`` and R installation.

    de_key: ``str``, optional, default: ``de_res:stat``
        Key name in the DE results to be used for gene signatures.
        Specify in format ``df_key:attr``, where ``df_key`` must exist in ``data.varm``, and ``attr`` must exist in columns of ``data.varm[df_key]``.

    minSize: ``int``, optional, default: ``5``
        Minimal size of a gene set to consider.

    maxSize: ``int``, optional, default: ``4000``
        Maximal size of a gene set to consider.

    n_jobs: ``int``, optional, default: ``1``
        Number of threads for parallel computation.
        If ``method='fgsea'``, set ``n_jobs=0`` by default to use one thread. If ``n_jobs>0``, set BPPARAM.

    seed: ``int``, optional, default: ``0``
        Random seed to make sure fGSEA results are reproducible.

    gsea_key: ``str``, optional, default: ``"gsea_out"``
        Key to use to store fGSEA results as a data frame.

    Returns
    -------
    ``None``

    Update ``data.uns``:
        ``data.uns[gsea_key]``: GSEA outputs sorted by ``padj``.

    Examples
    --------
    >>> pg.gsea(data, "de_res:stat", "hallmark", gsea_key='gsea_res')
    """
    if method == "blitzgsea":
        # center default is False, unless explicitly specified by users
        center = False
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
            center=center,
            **kwargs,
        )
    else:
        _run_fgsea(
            data=data,
            de_key=de_key,
            pathways=pathways,
            gsea_key=gsea_key,
            minSize=min_size,
            maxSize=max_size,
            nproc=n_jobs,
            seed=seed,
        )


def _split_de_key(de_key):
    assert ":" in de_key, f"Key '{de_key}' does not satisfy format 'df_key:attr'!"
    keys = de_key.split(":")
    df_key = keys[0]
    rank_key = ":".join(keys[1:])

    return df_key, rank_key


def _run_blitzgsea(
    data,
    de_key,
    pathways,
    gsea_key,
    min_size,
    max_size,
    n_jobs,
    seed,
    center,
    **kwargs,
) -> None:
    try:
        import blitzgsea as blitz
    except ModuleNotFoundError as e:
        import sys
        logger.error(f"{e}\nNeed blitzgsea! Try 'pip install blitzgsea'.")
        sys.exit(-1)

    df_key, rank_key = _split_de_key(de_key)
    assert df_key in data.varm, f"Key '{df_key}' not in data.varm! Wrong key name or run DE analysis beforehand!"
    assert rank_key in data.varm[df_key].dtype.names, f"Key '{rank_key}' not in DE result! Wrong key name specified!"

    rank_df = pd.DataFrame(data.varm[df_key], index=data.var_names)
    rank_df = rank_df.loc[~np.isnan(rank_df['padj'])].copy()
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
    minSize: int,
    maxSize: int,
    nproc: int,
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

    nproc = eff_n_jobs(nproc)
    if nproc == 1:
        nproc = 0  # Avoid BPPARAM to just use 1 core

    ro.r(f"set.seed({seed})")
    pwdict = load_signatures_from_file(predefined_pathways.get(pathways, pathways))
    pathways_r = ro.ListVector(pwdict)

    df_key, rank_key = _split_de_key(de_key)
    assert df_key in data.varm, f"Key '{df_key}' not in data.varm! Wrong key name or run DE analysis beforehand!"
    assert rank_key in data.varm[df_key].dtype.names, f"Key '{rank_key}' not in DE result! Wrong key name specified!"

    assert 'padj' in data.varm[df_key].dtype.names, f"No adjusted p-value exists in DE result!"
    qvals = data.varm[df_key]['padj']
    idx_select = np.where(~np.isnan(qvals))[0]  # Ignore genes with NaN adjusted p-values
    rank_vec = ro.FloatVector(data.varm[df_key][rank_key][idx_select])
    rank_vec.names = ro.StrVector(data.var_names[idx_select])
    res = fgsea.fgsea(pathways_r, rank_vec, minSize=minSize, maxSize=maxSize, nproc=nproc)
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
    res_df = res_df.set_index('Term')
    data.uns[gsea_key] = res_df


@timer(logger=logger)
def write_gsea_results_to_excel(
    data: Union[MultimodalData, UnimodalData],
    output_file: str,
    fgsea_key: Optional[str] = "fgsea_out",
    ndigits: Optional[int] = 3,
) -> None:
    """Write Gene Set Enrichment Analysis (GSEA) results generated by fgsea function into Excel workbook.

    Parameters
    ----------
    data: Union[``MultomodalData``, ``UnimodalData``]
        Single-cell or pseudo-bulk data.

    output_file: ``str``
        File name for the output.

    fgsea_key: ``str``, optinoal, default: ``"fgsea_out"``
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
