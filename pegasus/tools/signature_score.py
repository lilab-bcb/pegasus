import numpy as np
import pandas as pd

from typing import Dict, List
from anndata import AnnData

import logging

logger = logging.getLogger("pegasus")
from pegasus.utils import decorators as pg_deco


@pg_deco.TimeLogger()
def calc_signature_score(data: AnnData, signatures: Dict[str, List[str]], n_bins: int = 50) -> None:
    """Calculate signature / gene module score.

    This is an improved version of Livnat et al. 2018 Cell implementation.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    signatures: ``Dict[str, List[str]]``
        A dictionary containing multiple signature score calculation requests. Each score will be stored in data.obs field with key as the keyword. The value of the dict is a list of gene symbols.

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
    if "mean" not in data.var:
        data.var["mean"] = data.X.mean(axis = 0).A1

    if data.uns.get("sig_n_bins", 0) != n_bins:
        data.uns["sig_n_bins"] = n_bins
        mean_vec = data.var["mean"]
        bins = pd.qcut(mean_vec, n_bins)
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
        if idx.sum() < genes.size:
            omitted = ~genes.isin(data.var_names)
            logger.warning("For signature {}, genes {} are not in the data and thus omitted!".format(key, str(list(genes[omitted]))[1:-1]))
        if key in data.obs:
            logger.warning("Signature key {} exists in data.obs, the existing content will be overwritten!".format(key))
        data.obs[key] = (data.X[:, idx].mean(axis = 1).A1 - data.var.loc[idx, "mean"].mean()) - data.obsm["sig_background"][:, data.var["bins"].cat.codes[idx]].mean(axis = 1)
