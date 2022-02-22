import logging

logger = logging.getLogger(__name__)

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
