import numpy as np
import pandas as pd
import time
from natsort import natsorted

import multiprocessing
from sklearn.cluster import KMeans

from typing import List, Tuple
from pegasusio import UnimodalData, MultimodalData, read_input

import logging
logger = logging.getLogger(__name__)

def estimate_background_probs(hashing_data: UnimodalData, random_state: int = 0) -> None:
    """For cell-hashing data, estimate antibody background probability using KMeans algorithm.

    Parameters
    ----------
    hashing_data: ``UnimodalData``
        Annotated data matrix for antibody.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    Returns
    -------
    ``None``

    Update ``hashing_data``: Filtered cell barcodes with 0 counts

    Update ``hashing_data.uns``:
        * ``hashing_data.uns["background_probs"]``: estimated antibody background probability.

    Example
    -------
    >>> estimate_background_probs(hashing_data)
    """
    hashing_data.obs["counts"] = hashing_data.X.sum(axis=1).A1
    
    # Remove barcodes with 0 total counts
    idx = hashing_data.obs["counts"] == 0
    ncell_zero = idx.sum()
    if ncell_zero > 0:
        logger.warning(f"Detected {ncell_zero} cell barcodes with 0 hashtag counts, which are removed from the hashing data object.")
        hashing_data._inplace_subset_obs(~idx)        

    counts_log10 = np.log10(hashing_data.obs["counts"].values.reshape(-1, 1))
    kmeans = KMeans(n_clusters=2, random_state=random_state).fit(counts_log10)
    signal = 0 if kmeans.cluster_centers_[0] > kmeans.cluster_centers_[1] else 1
    hashing_data.obs["hto_type"] = "background"
    hashing_data.obs.loc[kmeans.labels_ == signal, "hto_type"] = "signal"

    idx = np.isin(hashing_data.obs["hto_type"], "background")
    pvec = hashing_data.X[idx,].sum(axis=0).A1
    back_probs = pvec / pvec.sum()

    idx = back_probs <= 0.0
    if idx.sum() > 0:
        logger.warning(f"Detected {idx.sum()} antibody barcodes {','.join(hashing_data.var_names[idx])} with 0 counts in the background! These barcodes are likely not in the experiment and thus removed.")
        hashing_data._inplace_subset_var(~idx)
        back_probs = back_probs[~idx]

    logger.info("Background probability distribution is estimated.")
    return back_probs


def estimate_probs(arr: List[float], pvec: List[float], alpha: float, alpha_noise: float, tol: float) -> List[float]:
    probs = np.zeros(pvec.size + 1)
    old_probs = np.zeros(pvec.size + 1)
    z = np.zeros(pvec.size + 1)
    noise = pvec.size

    # Estimate MLE without Generalized Dirichlet prior
    probs_mle = arr / arr.sum()
    probs[noise] = (probs_mle / pvec).min() + 0.01
    probs[:-1] = np.maximum(probs_mle - probs[noise] * pvec, 0.01)
    probs = probs / probs.sum()

    # EM algorithm
    i = 0
    eps = 1.0
    while eps > tol:
        i += 1
        old_probs[:] = probs[:]
        # E step
        z[:-1] = alpha - 1.0
        z[noise] = alpha_noise - 1.0
        for j in range(pvec.size):
            if arr[j] > 0:
                p = probs[j] / (probs[noise] * pvec[j] + probs[j])
                z[j] += arr[j] * p
                z[noise] += arr[j] * (1.0 - p)
        # M step
        idx = z > 0.0
        probs[idx] = z[idx] / z[idx].sum()
        probs[~idx] = 0.0
        eps = np.linalg.norm(probs - old_probs, ord=1)
        # print ("i = {}, eps = {:.2g}.".format(i, eps))

    return probs


def get_droplet_info(probs: List[float], sample_names: List[str]) -> Tuple[str, str]:
    ids = np.nonzero(probs >= 0.1)[0]
    ids = ids[np.argsort(probs[ids])[::-1]]
    return (
        "singlet" if ids.size == 1 else "doublet",
        ",".join([sample_names[i] for i in ids]),
    )


def calc_demux(raw_hto: UnimodalData, filt_hto: UnimodalData, nsample: int, min_signal: float, probs: str = "raw_probs") -> None:
    demux_type = np.full(raw_hto.shape[0], "unknown", dtype="object")
    assignments = np.full(raw_hto.shape[0], "", dtype="object")

    signals = filt_hto.obs["counts"].reindex(raw_hto.obs_names, fill_value=0.0).values * (
        1.0 - raw_hto.obsm[probs][:, nsample]
    )
    idx = signals >= min_signal

    tmp = raw_hto.obsm[probs][idx,]
    norm_probs = tmp[:, 0:nsample] / (1.0 - tmp[:, nsample])[:, None]

    values1 = []
    values2 = []
    for i in range(norm_probs.shape[0]):
        droplet_type, droplet_id = get_droplet_info(norm_probs[i,], filt_hto.var_names)
        values1.append(droplet_type)
        values2.append(droplet_id)

    demux_type[idx] = values1
    raw_hto.obs["demux_type"] = pd.Categorical(
        demux_type, categories=["singlet", "doublet", "unknown"]
    )
    assignments[idx] = values2
    raw_hto.obs["assignment"] = pd.Categorical(
        assignments, categories=natsorted(np.unique(assignments))
    )


def has_duplicate_names(names: List[str]) -> bool:
    for name in names:
        if name.find(".#~") >= 0:
            return True
    return False


def remove_suffix(assigns: List[str]) -> pd.Categorical:
    assigns = assigns.astype(str)
    results = []
    for value in assigns:
        fields = value.split(",")
        for i, item in enumerate(fields):
            pos = item.find(".#~")
            if pos >= 0:
                fields[i] = item[:pos]
        results.append(",".join(fields))
    results = np.array(results)
    return pd.Categorical(results, categories=natsorted(np.unique(results)))

def demultiplex(
    data: MultimodalData,
    selected: List[bool],
    min_signal: float = 10.0,
    alpha: float = 0.0,
    alpha_noise: float = 1.0,
    tol: float = 1e-6,
    n_threads: int = 1,
    **kwargs,
):
    """Demultiplexing cell/nucleus-hashing data, using the estimated antibody background probability calculated in ``demuxEM.estimate_background_probs``.

    Parameters
    ----------
    rna_data: ``UnimodalData``
        Data matrix for gene expression matrix.

    hashing_data: ``UnimodalData``
        Data matrix for HTO count matrix.

    min_signal: ``float``, optional, default: ``10.0``
        Any cell/nucleus with less than ``min_signal`` hashtags from the signal will be marked as ``unknown``.

    alpha: ``float``, optional, default: ``0.0``
        The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse.

    alpha_noise: ``float``, optional, default: ``1.0``
        The Dirichlet prior concenration parameter on the background noise.

    tol: ``float``, optional, default: ``1e-6``
        Threshold used for the EM convergence.

    n_threads: ``int``, optional, default: ``1``
        Number of threads to use. Must be a positive integer.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs["demux_type"]``: Demultiplexed types of the cells. Either ``singlet``, ``doublet``, or ``unknown``.
        * ``data.obs["assignment"]``: Assigned samples of origin for each cell barcode.
        * ``data.obs["assignment.dedup"]``: Only exist if one sample name can correspond to multiple feature barcodes. In this case, each feature barcode is assigned a unique sample name.

    Examples
    --------
    >>> demultiplex(hashing_data)
    """
    raw_hto = data.get_data(modality="hashing")
    prob_vector = estimate_background_probs(raw_hto, random_state=kwargs["random_state"])

    required_barcodes = pd.DataFrame(columns=['barcodekey', 'n_genes', 'n_counts'])
    required_barcodes.set_index('barcodekey')

    for i,required in enumerate(selected):
        if required == True:
            required_barcodes.append(raw_hto.obs.iloc[[i]])

    filt_hto = UnimodalData(required_barcodes,
                            raw_hto.var,
                            {"X":raw_hto.X},
                            None,
                            None,
                            None,
                            None,
                            None,
                            "X",
                            raw_hto.get_genome(),
                            raw_hto.get_modality(),
                            "filt_hto")

    filt_hto.uns["background_probs"] = prob_vector

    nsample = raw_hto.shape[1]
    assert (raw_hto.uns["background_probs"] <= 0.0).sum() == 0

    if nsample == 1:
        logger.warning("Detected only one barcode, no need to demultiplex!")
        raw_hto.obsm["raw_probs"] = np.zeros((filt_hto.shape[0], nsample + 1))
        raw_hto.obsm["raw_probs"][:, 0] = 1.0
        raw_hto.obsm["raw_probs"][:, 1] = 0.0
        raw_hto.obs["demux_type"] = "singlet"
        raw_hto.obs["assignment"] = raw_hto.var_names[0]
    else:
        if nsample == 2:
            logger.warning("Detected only two barcodes, demultiplexing accuracy might be affected!")
       
        ncalc = filt_hto.shape[0]
  
        hto_small = raw_hto.X.toarray()

        raw_hto.obsm["raw_probs"] = np.zeros((raw_hto.shape[0], nsample + 1))
        raw_hto.obsm["raw_probs"][:, nsample] = 1.0

        iter_array = [
            (hto_small[i,], raw_hto.uns["background_probs"], alpha, alpha_noise, tol)
            for i in range(ncalc)
        ]
        with multiprocessing.Pool(n_threads) as pool:
            raw_hto.obsm["raw_probs"] = pool.starmap(estimate_probs, iter_array)

        calc_demux(raw_hto, filt_hto, nsample, min_signal)

        if has_duplicate_names(raw_hto.var_names):
            raw_hto.obs["assignment.dedup"] = raw_hto.obs["assignment"]
            raw_hto.obs["assignment"] = remove_suffix(raw_hto.obs["assignment"].values)
    
    #data.update(hashing_data)
    data.update(filt_hto)
    data.select_data('filt_hto')
    logger.info("Demultiplexing is done.")