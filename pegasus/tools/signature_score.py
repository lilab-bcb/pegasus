import numpy as np
import pandas as pd

from typing import Dict, List, Union
from sklearn.cluster import KMeans

from pegasusio import UnimodalData, MultimodalData

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


import pkg_resources
predefined_signatures = dict(
    cell_cycle_human=pkg_resources.resource_filename("pegasus", "data_files/cell_cycle_human.gmt"),
    cell_cycle_mouse=pkg_resources.resource_filename("pegasus", "data_files/cell_cycle_mouse.gmt"),
    gender_human=pkg_resources.resource_filename("pegasus", "data_files/gender_human.gmt"),
    gender_mouse=pkg_resources.resource_filename("pegasus", "data_files/gender_mouse.gmt"),
    mitochondrial_genes_human=pkg_resources.resource_filename("pegasus", "data_files/mitochondrial_genes_human.gmt"),
    mitochondrial_genes_mouse=pkg_resources.resource_filename("pegasus", "data_files/mitochondrial_genes_mouse.gmt"),
    ribosomal_genes_human=pkg_resources.resource_filename("pegasus", "data_files/ribosomal_genes_human.gmt"),
    ribosomal_genes_mouse=pkg_resources.resource_filename("pegasus", "data_files/ribosomal_genes_mouse.gmt"),
)


def _check_and_calc_sig_background(data: UnimodalData, n_bins: int) -> bool:
    if "mean" not in data.var:
        data.var["mean"] = data.X.mean(axis = 0).A1

    if data.uns.get("sig_n_bins", 0) != n_bins:
        mean_vec = data.var["mean"]
        if mean_vec.size <= n_bins:
            logger.error(f"Number of bins {n_bins} is larger or equal to the total number of genes {mean_vec.size}! Please adjust n_bins and rerun this function!")
            return False
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

    return True


def _load_signatures_from_file(input_file: str) -> Dict[str, List[str]]:
    signatures = {}
    with open(input_file) as fin:
        for line in fin:
            items = line.strip().split('\t')
            signatures[items[0]] = list(set(items[2:]))
    logger.info(f"Loaded signatures from GMT file {input_file}.")
    return signatures


def _calc_sig_scores(data: UnimodalData, signatures: Dict[str, List[str]], show_omitted_genes: bool = False, skip_threshold: int = 2) -> None:
    for key, gene_list in signatures.items():
        genes = pd.Index(gene_list)
        idx = data.var_names.isin(genes)

        omit_string = ""
        nvalid = idx.sum()
        if nvalid < genes.size and show_omitted_genes:
            omitted = ~genes.isin(data.var_names)
            omit_string = f" Genes {str(list(genes[omitted]))[1:-1]} are not in the data and thus omitted."
        logger.info(f"Signature {key}: {nvalid} out of {genes.size} genes are used in signature score calculation.{omit_string}")

        if nvalid < skip_threshold:
            logger.warning(f"Signature {key} has less than {skip_threshold} genes kept and thus its score calculation is skipped!")
        else:
            if key in data.obs:
                logger.warning(f"Signature key {key} exists in data.obs, the existing content will be overwritten!")
            data.obs[key] = (data.X[:, idx].mean(axis = 1).A1 - data.var.loc[idx, "mean"].mean()) - data.obsm["sig_background"][:, data.var["bins"].cat.codes[idx]].mean(axis = 1)



@timer(logger=logger)
def calc_signature_score(data: MultimodalData, signatures: Union[Dict[str, List[str]], str], n_bins: int = 50, show_omitted_genes: bool = False, random_state: int = 0) -> None:
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

    if not _check_and_calc_sig_background(data, n_bins):
        return None

    if isinstance(signatures, str):
        sig_string = signatures
        if sig_string in predefined_signatures:
            signatures = _load_signatures_from_file(predefined_signatures[sig_string])
            if sig_string.startswith("cell_cycle"):
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes)
                kmeans = KMeans(n_clusters=3, random_state=random_state)
                kmeans.fit(data.obs[["G1/S", "G2/M"]].values)
                remap_code = np.zeros(3, dtype = np.int32)
                remap_code[1:] = np.argmax(kmeans.cluster_centers_, axis = 0)
                assert remap_code[1] != remap_code[2]
                remap_code[0] = 3 - remap_code[1:].sum()
                codes = remap_code[kmeans.labels_]
                data.obs["phase"] = pd.Categorical.from_codes(codes, categories = ["G0", "G1/S", "G2/M"])
            elif sig_string.startswith("gender"):
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
                data.obs["gender_score"] = data.obs["male_score"] - data.obs["female_score"]
                kmeans = KMeans(n_clusters=2, random_state=random_state)
                kmeans.fit(data.obs["gender_score"].values.reshape(-1, 1))
                codes = kmeans.labels_
                if kmeans.cluster_centers_[0, 0] > kmeans.cluster_centers_[1, 0]:
                    codes = 1 - codes # 0: female; 1: male
                data.obs["gender"] = pd.Categorical.from_codes(codes, categories = ["female", "male"])
            elif sig_string.startswith("mitochondrial_genes"):
                del signatures["mito_noncoding"]
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
            elif sig_string.startswith("ribosomal_genes"):
                del signatures["ribo_like"]
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
            else:
                assert False
        else:
            signatures = _load_signatures_from_file(sig_string)
            _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes)
    else:
        _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes)
