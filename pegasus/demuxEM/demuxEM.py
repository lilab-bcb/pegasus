import numpy as np
import pandas as pd
import time
from natsort import natsorted

import multiprocessing
from sklearn.cluster import KMeans

from typing import List
from anndata import AnnData

from pegasus.utils import decorators as pg_deco


def estimate_background_probs(adt: AnnData, random_state: int = 0):
    """For cell-hashing data, estimate antibody background probability using EM algorithm.

    Parameters
    ----------
    adt: ``anndata.AnnData``
        Annotated data matrix for antibody.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    Returns
    -------
    ``None``

    Update ``adt.uns``:
        * ``adt.uns["background_probs"]``: estimated antibody background probability.

    Example
    -------
    >>> pg.estimate_background_probs(adt)
    """
    adt.obs["counts"] = adt.X.sum(axis=1).A1 if adt.shape[1] > 1 else adt.X
    counts_log10 = np.log10(adt.obs["counts"].values.reshape(-1, 1))
    kmeans = KMeans(n_clusters=2, random_state=random_state).fit(counts_log10)
    signal = 0 if kmeans.cluster_centers_[0] > kmeans.cluster_centers_[1] else 1
    adt.obs["hto_type"] = "background"
    adt.obs.loc[kmeans.labels_ == signal, "hto_type"] = "signal"

    idx = np.isin(adt.obs["hto_type"], "background")
    pvec = (
        adt.X[idx,].sum(axis=0).A1 if adt.shape[1] > 1 else np.array(adt.X[idx,].sum())
    )
    pvec /= pvec.sum()

    adt.uns["background_probs"] = pvec


def estimate_probs(arr, pvec, alpha, alpha_noise, tol):
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


def get_droplet_info(probs, sample_names):
    ids = np.nonzero(probs >= 0.1)[0]
    ids = ids[np.argsort(probs[ids])[::-1]]
    return (
        "singlet" if ids.size == 1 else "doublet",
        ",".join([sample_names[i] for i in ids]),
    )


def calc_demux(data, adt, nsample, min_signal, probs="raw_probs"):
    demux_type = np.full(data.shape[0], "unknown", dtype="object")
    assignments = np.full(data.shape[0], "", dtype="object")

    signals = adt.obs["counts"].reindex(data.obs_names, fill_value=0.0).values * (
        1.0 - data.obsm[probs][:, nsample]
    )
    idx = signals >= min_signal

    tmp = data.obsm[probs][idx,]
    norm_probs = tmp[:, 0:nsample] / (1.0 - tmp[:, nsample])[:, None]

    values1 = []
    values2 = []
    for i in range(norm_probs.shape[0]):
        droplet_type, droplet_id = get_droplet_info(norm_probs[i,], adt.var_names)
        values1.append(droplet_type)
        values2.append(droplet_id)

    demux_type[idx] = values1
    data.obs["demux_type"] = pd.Categorical(
        demux_type, categories=["singlet", "doublet", "unknown"]
    )
    assignments[idx] = values2
    data.obs["assignment"] = pd.Categorical(
        assignments, categories=natsorted(np.unique(assignments))
    )

def has_duplicate_names(names: List[str]) -> bool:
    for name in names:
        if name.find(".#~") >= 0:
            return True
    return False

def remove_suffix(assigns: List[str]) -> List[str]:
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


@pg_deco.TimeLogger()
def demultiplex(
    data: AnnData,
    adt: AnnData,
    min_signal: float = 10.0,
    alpha: float = 0.0,
    alpha_noise: float = 1.0,
    tol: float = 1e-6,
    n_threads: int = 1,
):
    """Demultiplexing cell-hashing data, using the estimated antibody background probability calculated in ``pg.estimate_background_probs``.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix for gene expression matrix.

    adt: ``anndata.AnnData``
        Annotated data matrix for antibody count matrix.

    min_signal: ``float``, optional, default: ``10.0``
        Any cell/nucleus with less than ``min_signal`` hashtags from the signal will be marked as ``unknown``.

    alpha: ``float``, optional, default: ``0.0``
        The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse.

    alpha_noise: ``float``, optional, default: ``1.0``

    tol: ``float``, optional, default: ``1e-6``
        Tolerance threshold when judging equivalence of two floating point values.

    n_threads: ``int``, optional, default: ``1``
        Number of threads to use. Must be a positive integer.

    Returns
    -------
    ``None``

    Update ``data.obs``:
        * ``data.obs["demux_type"]``: Demultiplexed types of the cells. Either ``singlet``, ``doublet``, or ``unknown``.
        * ``data.obs["assignment"]``: Assigned samples of origin for each cell barcode.
        * ``data.obs["assignment.dedup"]``: Only exist if one sample name can correspond to multiple feature barcodes. In this array, each feature barcode is assigned a unique sample name.

    Examples
    --------
    >>> pg.demultiplex(adata, adt)
    """

    nsample = adt.shape[1]
    data.uns["background_probs"] = adt.uns["background_probs"]

    idx_df = data.obs_names.isin(adt.obs_names)
    adt.obs["rna_type"] = "background"
    adt.obs.loc[data.obs_names[idx_df], "rna_type"] = "signal"

    if nsample == 1:
        print("Warning: detect only one barcode, no need to demultiplex!")
        data.obsm["raw_probs"] = np.zeros((data.shape[0], nsample + 1))
        data.obsm["raw_probs"][:, 0] = 1.0
        data.obsm["raw_probs"][:, 1] = 0.0
        data.obs["demux_type"] = "singlet"
        data.obs["assignment"] = adt.var_names[0]
    else:
        if nsample == 2:
            print(
                "Warning: detect only two barcodes, demultiplexing accuracy might be affected!"
            )

        ncalc = idx_df.sum()
        if ncalc < data.shape[0]:
            nzero = data.shape[0] - ncalc
            print(
                "Warning: {} cells do not have ADTs, percentage = {:.2f}%.".format(
                    nzero, nzero * 100.0 / data.shape[0]
                )
            )
        adt_small = adt[data.obs_names[idx_df],].X.toarray()

        data.obsm["raw_probs"] = np.zeros((data.shape[0], nsample + 1))
        data.obsm["raw_probs"][:, nsample] = 1.0

        iter_array = [
            (adt_small[i,], adt.uns["background_probs"], alpha, alpha_noise, tol)
            for i in range(ncalc)
        ]
        with multiprocessing.Pool(n_threads) as pool:
            data.obsm["raw_probs"][idx_df, :] = pool.starmap(estimate_probs, iter_array)

        calc_demux(data, adt, nsample, min_signal)

        if has_duplicate_names(adt.var_names):
            data.obs["assignment.dedup"] = data.obs["assignment"]
            data.obs["assignment"] = remove_suffix(data.obs["assignment"].values)

