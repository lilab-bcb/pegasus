import numpy as np
import pandas as pd

from typing import Dict, List, Union
from sklearn.cluster import KMeans

import anndata
from pegasusio import UnimodalData, MultimodalData
from pegasus.tools import calc_mean, calc_sig_background

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
    apoptosis_human=pkg_resources.resource_filename("pegasus", "data_files/apoptosis_human.gmt"),
    apoptosis_mouse=pkg_resources.resource_filename("pegasus", "data_files/apoptosis_mouse.gmt"),
)


def _check_and_calc_sig_background(data: UnimodalData, n_bins: int) -> bool:
    if "mean" not in data.var:
        data.var["mean"] = calc_mean(data.X, axis = 0)

    if data.uns.get("sig_n_bins", 0) != n_bins:
        mean_vec = data.var["mean"].values
        if mean_vec.size <= n_bins:
            logger.error(f"Number of bins {n_bins} is larger or equal to the total number of genes {mean_vec.size}! Please adjust n_bins and rerun this function!")
            return False
        data.uns["sig_n_bins"] = n_bins
        try:
            bins = pd.qcut(mean_vec, n_bins)
        except ValueError:
            logger.warning("Detected and dropped duplicate bin edges!")
            bins = pd.qcut(mean_vec, n_bins, duplicates = "drop")
        if bins.value_counts().min() == 1:
            logger.warning("Detected bins with only 1 gene!")
        bins.categories = bins.categories.astype(str)
        data.var["bins"] = bins

        # calculate background expectations
        data.obsm["sig_bkg_mean"], data.obsm["sig_bkg_std"] = calc_sig_background(data.X, bins, mean_vec)

    return True


def _load_signatures_from_file(input_file: str) -> Dict[str, List[str]]:
    signatures = {}
    with open(input_file) as fin:
        for line in fin:
            items = line.strip().split('\t')
            signatures[items[0]] = list(set(items[2:]))
    logger.info(f"Loaded signatures from GMT file {input_file}.")
    return signatures


def _calc_sig_scores(data: UnimodalData, signatures: Dict[str, List[str]], show_omitted_genes: bool = False, skip_threshold: int = 1) -> None:
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

            data.obs[key] = ((data.X[:, idx].toarray() - data.var.loc[idx, "mean"].values - data.obsm["sig_bkg_mean"][:, data.var["bins"].cat.codes[idx]]) / data.obsm["sig_bkg_std"][:,data.var["bins"].cat.codes[idx]]).mean(axis = 1).astype(np.float32)
            data.register_attr(key, "signature")


def calculate_z_score(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    n_bins: int = 50,
) -> np.array:
    """Calculate the standardized z scores of the count matrix.

    Parameters
    -----------
    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData`` object.
        Single cell expression data.
    n_bins: ``int``, optional, default: ``50``
        Number of bins on expression levels for grouping genes.

    Returns
    -------
    numpy.array
        A 2D numpy array of shape ``(n_cells, n_features)``, which represents the standardized z-score expression matrix.

    Examples
    ---------
    >>> pg.calculate_z_score(data)
    >>> pg.calculate_z_score(data, n_bins=100)
    """
    if isinstance(data, MultimodalData):
        data = data._unidata

    if not _check_and_calc_sig_background(data, n_bins):
        return None

    z_score_mat = (data.X.toarray().astype(np.float32) - data.var["mean"].values.astype(np.float32) - data.obsm["sig_bkg_mean"][:, data.var["bins"].cat.codes].astype(np.float32)) / data.obsm["sig_bkg_std"][:, data.var["bins"].cat.codes].astype(np.float32)

    return z_score_mat


@timer(logger=logger)
def calc_signature_score(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    signatures: Union[Dict[str, List[str]], str],
    n_bins: int = 50,
    show_omitted_genes: bool = False,
    random_state: int = 0
) -> None:
    """Calculate signature / gene module score. [Li20-1]_

    This is an improved version of implementation in [Jerby-Arnon18]_.

    Parameters
    ----------
    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData`` object.
        Single cell expression data.

    signatures: ``Dict[str, List[str]]`` or ``str``
        This argument accepts either a dictionary or a string.
        If ``signatures`` is a dictionary, it can contain multiple signature score calculation requests. Each key in the dictionary represents a separate signature score calculation and its corresponding value contains a list of gene symbols. Each score will be stored in data.obs field with key as the keyword.
        If ``signatures`` is a string, it should refer to a Gene Matrix Transposed (GMT)-formatted file. Pegasus will load signatures from the GMT file.

        Pegasus also provide 5 default signature panels for each of human and mouse. They are ``cell_cycle_human``, ``gender_human``, ``mitochondrial_genes_human``, ``ribosomal_genes_human`` and ``apoptosis_human`` for human; ``cell_cycle_mouse``, ``gender_mouse``, ``mitochondrial_genes_mouse``, ``ribosomal_genes_mouse`` and ``apoptosis_mouse`` for mouse.
            * ``cell_cycle_human`` contains two cell-cycle signatures, ``G1/S`` and ``G2/M``, obtained from Tirosh et al. 2016. We also updated gene symbols according to Seurat's ``cc.genes.updated.2019`` vector. We additionally calculate signature scores ``cycling`` and ``cycle_diff``, which are ``max{G2/M, G1/S}`` and ``G2/M`` - ``G1/S`` respectively. We provide predicted cell cycle phases in ``data.obs['predicted_phase']`` in case it is useful. ``predicted_phase`` is predicted based on ``G1/S`` and ``G2/M`` scores. First, we identify ``G0`` cells. We apply KMeans algorithm to obtain 2 clusters based on the ``cycling`` signature. ``G0`` cells are from the cluster with smallest mean value. For each cell from the other cluster, if ``G1/S`` > ``G2/M``, it is a ``G1/S`` cell, otherwise it is a ``G2/M`` cell.
            * ``gender_human`` contains two gender-specific signatures, ``female_score`` and ``male_score``. Genes were selected based on DE analysis between genders based on 8 channels of bone marrow data from HCA Census of Immune Cells and the brain nuclei data from Gaublomme and Li et al, 2019, Nature Communications. After calculation, three signature scores will be calculated: ``female_score``, ``male_score`` and ``gender_score``. ``female_score`` and ``male_score`` are calculated based on female and male signatures respectively and a larger score represent a higher likelihood of that gender. ``gender_score`` is calculated as ``male_score`` - ``female_score``. A large positive score likely represents male and a large negative score likely represents female. Pegasus also provides predicted gender for each cell based on ``gender_score``, which is stored in ``data.obs['predicted_gender']``. To predict genders, we apply the KMeans algorithm to the ``gender_score`` and ask for 3 clusters. The clusters with a minimum and maximum clauster centers are predicted as ``female`` and ``male`` respectively and the cluster in the middle is predicted as ``uncertain``. Note that this approach is conservative and it is likely that users can predict genders based on ``gender_score`` for cells in the ``uncertain`` cluster with a reasonable accuracy.
            * ``mitochondrial_genes_human`` contains two signatures, ``mito_genes`` and ``mito_ribo``. ``mito_genes`` contains 13 mitocondrial genes from chrM and ``mito_ribo`` contains mitocondrial ribosomal genes that are not from chrM. Note that ``mito_genes`` correlates well with percent of mitocondrial UMIs and ``mito_ribo`` does not.
            * ``ribosomal_genes_human`` contains one signature, ``ribo_genes``, which includes ribosomal genes from both large and small units.
            * ``apoptosis_human`` contains one signature, ``apoptosis``, which includes apoptosis-related genes from the KEGG pathway.
            * ``cell_cycle_mouse``, ``gender_mouse``, ``mitochondrial_genes_mouse``, ``ribosomal_genes_mouse`` and ``apoptosis_mouse`` are the corresponding signatures for mouse. Gene symbols are directly translated from human genes.

    n_bins: ``int``, optional, default: 50
        Number of bins on expression levels for grouping genes.

    show_omitted_genes: ``bool``, optional, default False
        Signature genes that are not expressed in the data will be omitted. By default, pegasus does not report which genes are omitted. If this option is turned on, report omitted genes.

    random_state: ``int``, optional, default: 0
        Random state used by KMeans if signature == ``gender_human`` or ``gender_mouse``.

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
    >>> pg.calc_signature_score(data, "cell_cycle_human")
    """
    if isinstance(data, MultimodalData):
        data = data._unidata

    if not _check_and_calc_sig_background(data, n_bins):
        return None

    if isinstance(signatures, str):
        sig_string = signatures
        if sig_string in predefined_signatures:
            signatures = _load_signatures_from_file(predefined_signatures[sig_string])
            from threadpoolctl import threadpool_limits

            if sig_string.startswith("cell_cycle"):
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes)
                data.obs["cycle_diff"] = data.obs["G2/M"] - data.obs["G1/S"]

                values = data.obs[["G1/S", "G2/M"]].values
                maxvalues = values.max(axis = 1)
                data.obs["cycling"] = maxvalues

                kmeans = KMeans(n_clusters=2, random_state=random_state)
                with threadpool_limits(limits = 1):
                    kmeans.fit(maxvalues.reshape(-1, 1))
                cycle_idx = kmeans.labels_ == np.argmax(kmeans.cluster_centers_[:,0])

                codes = np.zeros(data.shape[0], dtype = np.int32)
                codes[cycle_idx & (values[:, 0] == maxvalues)] = 1
                codes[cycle_idx & (values[:, 1] == maxvalues)] = 2

                data.obs["predicted_phase"] = pd.Categorical.from_codes(codes, categories = ["G0", "G1/S", "G2/M"])
            elif sig_string.startswith("gender"):
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
                data.obs["gender_score"] = data.obs["male_score"] - data.obs["female_score"]

                kmeans = KMeans(n_clusters=3, random_state=random_state)
                with threadpool_limits(limits = 1):
                    kmeans.fit(data.obs["gender_score"].values.reshape(-1, 1))
                reorg_dict = {y:x for x, y in enumerate(np.argsort(kmeans.cluster_centers_[:,0]))}
                codes = list(map(lambda x: reorg_dict[x], kmeans.labels_))

                data.obs["predicted_gender"] = pd.Categorical.from_codes(codes, categories = ["female", "uncertain", "male"])
            elif sig_string.startswith("mitochondrial_genes"):
                del signatures["mito_noncoding"]
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
            elif sig_string.startswith("ribosomal_genes"):
                del signatures["ribo_like"]
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
            elif sig_string.startswith("apoptosis"):
                _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes, skip_threshold = 1)
            else:
                assert False
        else:
            signatures = _load_signatures_from_file(sig_string)
            _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes)
    else:
        _calc_sig_scores(data, signatures, show_omitted_genes = show_omitted_genes)
