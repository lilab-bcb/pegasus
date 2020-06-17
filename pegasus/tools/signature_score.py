import numpy as np
import pandas as pd

from typing import Dict, List, Union
from pegasusio import UnimodalData, MultimodalData

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer



@timer(logger=logger)
def calc_signature_score(data: MultimodalData, signatures: Union[Dict[str, List[str]], str], n_bins: int = 50) -> None:
    """Calculate signature / gene module score.

    This is an improved version of Livnat et al. 2018 Cell implementation.

    Parameters
    ----------
    data: ``anndata.AnnData``
        MultimodalData with the current selected UnimodalData used.

    signatures: ``Dict[str, List[str]]`` or ``str``
        A dictionary containing multiple signature score calculation requests. Each score will be stored in data.obs field with key as the keyword. The value of the dict is a list of gene symbols. If <signatures> is a string, load signatures from the corresponding Gene Matrix Transposed (GMT)-formatted file. 

    n_bins: ``int``, optional, default: 50

    Returns
    -------
    ``None``.

    Update ``data.obs``:

        * ``data.obs["key"]``: signature / gene module score for signature "key"

    Update ``data.var``:

        * ``data.var["mean"]``: Mean expression of each gene across all cells. Only updated if "mean" does not exist in data.var.

        * ``data.var["bins"]``: Bin category for each gene. Only updated if data.uns["sig_n_bins"] is updated.

    Update ``data.obsm``:

        * ``data.obsm["sig_background"]``: Expected signature score for each bin category. Only updated if data.uns["sig_n_bins"] is updated.

    Update ``data.uns``:

        * ``data.uns["sig_n_bins"]``: Number of bins to partition genes into. Only updated if "sig_n_bins" does not exist or the recorded number of bins does not match n_bins.

    Examples
    --------
    >>> pg.calc_signature_score(data, {"T_cell_sig": ["CD3D", "CD3E", "CD3G", "TRAC"]})
    """
    if isinstance(data, MultimodalData):
        data = data._unidata

    if isinstance(signatures, str):
        input_file = signatures
        signatures = {}
        with open(input_file) as fin:
            for line in fin:
                items = line.strip().split('\t')
                signatures[items[0]] = list(set(items[2:]))
        logger.info(f"Loaded signatures from GMT file {input_file}.")


    if "mean" not in data.var:
        data.var["mean"] = data.X.mean(axis = 0).A1

    if data.uns.get("sig_n_bins", 0) != n_bins:
        mean_vec = data.var["mean"]
        if mean_vec.size <= n_bins:
            logger.error(f"Number of bins {n_bins} is larger or equal to the total number of genes {mean_vec.size}! Please adjust n_bins and rerun this function!")
            return None
        data.uns["sig_n_bins"] = n_bins
        try:
            bins = pd.qcut(mean_vec, n_bins)
        except ValueError:
            logger.warning("Detected and dropped duplicate bin edges!")
            bins = pd.qcut(mean_vec, n_bins, duplicates = "drop")
            n_bins = bins.cat.categories.size
        if bins.value_counts().min() == 1:
            logger.warning("Detected bins with only 1 gene!")
        bins.cat.categories = bins.cat.categories.astype(str)
        data.var["bins"] = bins
        # calculate background expectations
        sig_background = np.zeros((data.shape[0], n_bins))
        for code in range(n_bins):
            idx = (bins.cat.codes == code).values
            base = mean_vec[idx].mean()
            sig_background[:, code] = data.X[:, idx].mean(axis = 1).A1 - base
        data.obsm["sig_background"] = sig_background

    for key, gene_list in signatures.items():
        genes = pd.Index(gene_list)
        idx = data.var_names.isin(genes)
        nvalid = idx.sum()
        if nvalid < genes.size:
            omitted = ~genes.isin(data.var_names)
            logger.warning(f"For signature {key}, genes {str(list(genes[omitted]))[1:-1]} are not in the data and thus omitted!")
        if nvalid < 2:
            logger.warning(f"Signature {key} has less than 2 genes in the data and thus we skip its score calculation!")
        else:
            if key in data.obs:
                logger.warning("Signature key {} exists in data.obs, the existing content will be overwritten!".format(key))
            data.obs[key] = (data.X[:, idx].mean(axis = 1).A1 - data.var.loc[idx, "mean"].mean()) - data.obsm["sig_background"][:, data.var["bins"].cat.codes[idx]].mean(axis = 1)
