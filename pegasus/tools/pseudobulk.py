import numpy as np
import pandas as pd
import logging
logger = logging.getLogger(__name__)

from pegasusio import MultimodalData, UnimodalData, timer
from pegasus.tools import eff_n_jobs
from pandas.api.types import is_numeric_dtype
from typing import Union, Optional, List, Tuple



def set_bulk_value(col):
    if is_numeric_dtype(col):
        return col.mean()
    else:
        # Categorical
        return col.value_counts().idxmax()


def get_pseudobulk_count(X, df, attr, bulk_list):
    M = []
    for bulk in bulk_list:
        df_bulk = df.loc[df[attr] == bulk]
        bulk_vec = np.sum(X[df_bulk.index, :], axis=0).A1
        M.append(bulk_vec)

    return np.array(M, dtype=np.int32)


@timer(logger=logger)
def pseudobulk(
    data: Union[MultimodalData, UnimodalData],
    groupby: str,
    attrs: Optional[Union[List[str], str]] = None,
    mat_key: Optional[str] = None,
    condition: Optional[str] = None,
) -> MultimodalData:
    """Generate Pseudo-bulk count matrices.

    Parameters
    -----------
    data: ``MultimodalData`` or ``UnimodalData`` object
        Annotated data matrix with rows for cells and columns for genes.

    groupby: ``str``
        Specify the cell attribute used for aggregating pseudo-bulk data.
        Key must exist in ``data.obs``.

    attrs: ``str`` or ``List[str]``, optional, default: ``None``
        Specify additional cell attributes to remain in the pseudo bulk data.
        If set, all attributes' keys must exist in ``data.obs``.
        Notice that for a categorical attribute, each pseudo-bulk's value is the one of highest frequency among its cells,
        and for a numeric attribute, each pseudo-bulk's value is the mean among its cells.

    mat_key: ``str``, optional, default: ``None``
        Specify the single-cell count matrix used for aggregating pseudo-bulk counts:
        If specified, use the count matrix with key ``mat_key`` from matrices of ``data``; otherwise, first look for key ``counts``, then for ``raw.X`` if not existing.

    condition: ``str``, optional, default: ``None``
        If set, additionally generate pseudo-bulk matrices per condition specified in ``data.obs[condition]``.

    Returns
    -------
    A MultimodalData object ``mdata`` containing pseudo-bulk information:
        * It has the following count matrices:

          * ``X``: The pseudo-bulk count matrix over all cells.
          * If ``condition`` is set, add additional pseudo-bulk count matrices of cells restricted to each condition, respectively
        * ``mdata.obs``: It contains pseudo-bulk attributes aggregated from the corresponding single-cell attributes.
        * ``mdata.var``: Gene names and Ensembl IDs are maintained.

    Examples
    --------
    >>> pg.pseudobulk(data, groupby="Channel")
    """
    if mat_key is None:
        if "counts" in data._unidata.matrices:
            mat_key = "counts"
        elif "raw.X" in data._unidata.matrices:
            mat_key = "raw.X"
        else:
            import sys
            logger.error("No matrix with default key found in data! Please specify an explicit matrix key!")
            sys.exit(-1)
    X = data.get_matrix(mat_key)

    if groupby not in data.obs.columns:
        import sys
        logger.error(f"Sample key '{groupby}' must exist in data.obs!")
        sys.exit(-1)

    sample_vec = (
        data.obs[groupby]
        if isinstance(data.obs[groupby].dtype, pd.CategoricalDtype)
        else data.obs[groupby].astype("category")
    )
    bulk_list = sample_vec.cat.categories

    df_barcode = data.obs.reset_index()

    mat_dict = {"counts": get_pseudobulk_count(X, df_barcode, groupby, bulk_list)}

    # Generate pseudo-bulk attributes if specified
    bulk_attr_list = []

    if attrs is not None:
        if isinstance(attrs, str):
            attrs = [attrs]
        for attr in attrs:
            if attr not in data.obs.columns:
                import sys
                logger.error(f"Cell attribute key '{attr}' must exist in data.obs!")
                sys.exit(-1)

    for bulk in bulk_list:
        df_bulk = df_barcode.loc[df_barcode[groupby] == bulk]
        if attrs is not None:
            bulk_attr = df_bulk[attrs].apply(set_bulk_value, axis=0)
            bulk_attr["barcodekey"] = bulk
        else:
            bulk_attr = pd.Series({"barcodekey": bulk})
        bulk_attr_list.append(bulk_attr)

    df_pseudobulk = pd.DataFrame(bulk_attr_list)
    for col in df_pseudobulk.columns:
        if col == 'barcodekey':
            continue
        if isinstance(df_barcode[col].dtype, pd.CategoricalDtype):
            c_old = df_barcode[col].cat
            cat_dtype = pd.CategoricalDtype(categories=c_old.categories, ordered=c_old.ordered)
            c_new = pd.Series(df_pseudobulk[col], dtype=cat_dtype).cat
            df_pseudobulk[col] = c_new.remove_unused_categories()

    df_feature = pd.DataFrame(index=data.var_names)
    if "featureid" in data.var.columns:
        df_feature["featureid"] = data.var["featureid"]

    if condition is not None:
        if condition not in data.obs.columns:
            import sys
            logger.error(f"Condition key '{attr}' must exist in data.obs!")
            sys.exit(-1)

        cluster_list = data.obs[condition].astype("category").cat.categories
        for cls in cluster_list:
            mat_dict[f"{condition}_{cls}.X"] = get_pseudobulk_count(
                X, df_barcode.loc[df_barcode[condition] == cls], groupby, bulk_list
            )

    udata = UnimodalData(
        barcode_metadata=df_pseudobulk,
        feature_metadata=df_feature,
        matrices=mat_dict,
        genome=data.get_genome(),
        modality="pseudobulk",
        cur_matrix="counts",
    )

    return MultimodalData(udata)


@timer(logger=logger)
def deseq2(
    pseudobulk: MultimodalData,
    design: str,
    contrasts: Union[Tuple[str, str, str], List[Tuple[str, str, str]]],
    backend: str = "pydeseq2",
    de_key: Union[str, List[str]] = "deseq2",
    alpha: float = 0.05,
    compute_all: bool = False,
    verbose: bool = True,
    n_jobs: int = -1,
) -> None:
    """Perform Differential Expression (DE) Analysis using DESeq2 on pseduobulk data.

    Parameters
    ----------
    pseudobulk: ``MultimodalData`` or ``UnimodalData``
        Pseudobulk data with rows for samples/pseudobulks and columns for genes. It may contain multiple count matrices of the same shape with different keys

    design: ``str`` or ``List[str]``
        For ``pydeseq2`` backend, specify either a factor or a list of factors to be used as design variables.They must be all in ``pseudobulk.obs``.
        For ``deseq2`` backend, specify the design formula that will be passed to DESeq2. E.g. ``~group+condition`` or ``~genotype+treatment+genotype:treatment``.

    contrasts: ``Tuple[str, str, str]`` or ``List[Tuple[str, str, str]]``
        A tuple of three elements passing to DESeq2: a factor in design formula, a level in the factor as the test level (numeritor of fold change), and a level as the reference level (denominator of fold change).
        It also accept multiple contrasts as a list. In this way, ``de_key`` must be a list of strings, and the DE result of each contrast will then be stored in ``data.varm`` with the corresponding key.

    backend: ``str``, optional, default: ``pydeseq2``
        Specify which package to use as the backend for pseudobulk DE analysis.
        By default, use ``PyDESeq2`` which is a pure Python implementation of DESeq2 method.
        Alternatively, if specifying ``deseq2``, then use R package DESeq2, which requires ``rpy2`` package and R installation.

    de_key: ``str`` or ``List[str]``, optional, default: ``"deseq2"``
        Key name of DE analysis results stored. For count matrix with name ``condition.X``, stored key will be ``condition.de_key``.
        Provide a list of keys if ``contrasts`` is a list.

    alpha: ``float``, optional, default: ``0.05``
        The significance cutoff (between 0 and 1) used for optimizing the independent filtering to calculate the adjusted p-values (FDR).

    compute_all: ``bool``, optional, default: ``False``
        If performing DE analysis on all count matrices. By default (``compute_all=False``), only apply DE analysis to the default count matrix ``counts``.

    verbose: ``bool``, optional, default: ``True``
        If showing DESeq2 status updates during fit. Only works when ``backend="pydeseq2"``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to use. If ``-1``, use all physical CPU cores. This only works when ``backend="pydeseq2"`.

    Returns
    -------
    ``None``

    Update ``pseudobulk.varm``:
        ``pseudobulk.varm[de_key]``: DE analysis result for pseudo-bulk count matrix.
        (Optional) ``pseudobulk.varm[condition.de_key]``: If ``compute_all=True``, DE results for each condition-specific pseudo-bulk count matrices.

    Examples
    --------
    >>> pg.deseq2(pseudobulk, 'gender', ('gender', 'female', 'male'))
    >>> pg.deseq2(pseudobulk, '~gender', ('gender', 'female', 'male'), backend="deseq2")
    >>> pg.deseq2(pseudobulk, 'treatment', [('treatment', 'A', 'B'), ('treatment', 'A', 'C')], de_key=['deseq2_A_B', 'deseq2_A_C'])
    """
    mat_keys = ['counts'] if not compute_all else pseudobulk.list_keys()
    for mat_key in mat_keys:
        if backend == "pydeseq2":
            _run_pydeseq2(pseudobulk=pseudobulk, mat_key=mat_key, design_factors=design, contrasts=contrasts, de_key=de_key, alpha=alpha, n_jobs=n_jobs, verbose=verbose)
        else:
            _run_rdeseq2(pseudobulk=pseudobulk, mat_key=mat_key, design=design, contrasts=contrasts, de_key=de_key, alpha=alpha)


def _run_pydeseq2(
    pseudobulk: MultimodalData,
    mat_key: str,
    design_factors: Union[str, List[str]],
    contrasts: Union[Tuple[str, str, str], List[Tuple[str, str, str]]],
    de_key: Union[str, List[str]],
    alpha: float,
    n_jobs: int,
    verbose: bool,
) -> None:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
        from pydeseq2.ds import DeseqStats
    except ModuleNotFoundError as e:
        import sys
        logger.error(f"{e}\nNeed pydeseq2! Try 'pip install pydeseq2'.")
        sys.exit(-1)

    if isinstance(design_factors, str):
        if design_factors not in pseudobulk.obs.columns:
            import sys
            logger.error(f"The design factor {design_factors} does not exist in data.obs!")
            sys.exit(-1)
    else:
        for factor in design_factors:
            if factor not in pseudobulk.obs.columns:
                import sys
                logger.error(f"The design factor {factor} does not exist in data.obs!")
                sys.exit(-1)

    counts_df = pd.DataFrame(pseudobulk.get_matrix(mat_key), index=pseudobulk.obs_names, columns=pseudobulk.var_names)
    metadata = pseudobulk.obs

    n_cpus = eff_n_jobs(n_jobs)
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=design_factors,
        inference=inference,
        quiet=not verbose,
    )
    dds.deseq2()

    if isinstance(contrasts, Tuple):
        contrasts = [contrasts]
    if isinstance(de_key, str):
        de_key = [de_key]

    for i, contrast in enumerate(contrasts):
        stat_res = DeseqStats(
            dds,
            contrast=contrast,
            alpha=alpha,
            inference=inference,
            quiet=not verbose,
        )
        stat_res.summary()
        res_key = de_key[i] if mat_key == "counts" else mat_key.removesuffix(".X") + "." + de_key[i]
        res_df = stat_res.results_df
        pseudobulk.varm[res_key] = res_df.to_records(index=False)


def _run_rdeseq2(
    pseudobulk: MultimodalData,
    mat_key: str,
    design: str,
    contrasts: Union[Tuple[str, str, str], List[Tuple[str, str, str]]],
    de_key: Union[str, List[str]],
    alpha: float,
) -> None:
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, numpy2ri, Formula
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter
    except ModuleNotFoundError as e:
        import sys
        logger.error(f"{e}\nNeed rpy2! Try 'pip install rpy2'.")
        sys.exit(-1)

    try:
        deseq2 = importr('DESeq2')
    except ModuleNotFoundError:
        import sys
        text = """Please install DESeq2 in order to run this function.\n
                To install this package, start R and enter:\n
                if (!require("BiocManager", quietly = TRUE))
                    install.packages("BiocManager")
                BiocManager::install("DESeq2")"""

        logger.error(text)
        sys.exit(-1)

    if not design.strip().startswith("~"):
        import sys
        logger.error(f"Design '{design}' is not a valid R formula! Valid examples: '~var', '~var1+var2', '~var1+var2+var1:var2'.")
        sys.exit(-1)

    import math
    to_dataframe = ro.r('function(x) data.frame(x)')

    with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):
        dds = deseq2.DESeqDataSetFromMatrix(countData = pseudobulk.get_matrix(mat_key).T, colData = pseudobulk.obs, design = Formula(design))

    dds = deseq2.DESeq(dds)

    if isinstance(contrasts, Tuple):
        contrasts = [contrasts]
    if isinstance(de_key, str):
        de_key = [de_key]

    for i, contrast in enumerate(contrasts):
        res= deseq2.results(dds, contrast=ro.StrVector(contrast), alpha=alpha)
        with localconverter(ro.default_converter + pandas2ri.converter):
            res_df = ro.conversion.rpy2py(to_dataframe(res))

        res_key = de_key[i] if mat_key == "counts" else mat_key.removesuffix(".X") + "." + de_key[i]
        pseudobulk.varm[res_key] = res_df.to_records(index=False)
