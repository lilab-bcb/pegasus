import logging

logger = logging.getLogger(__name__)

import os
import itertools
import scipy.sparse
import numpy as np
import pandas as pd

from natsort import natsorted
from scipy.sparse import csr_matrix, issparse
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
from pegasus.tools import predefined_gene_orders, process_mat_key, log_norm, eff_n_jobs
from pegasusio import MultimodalData, UnimodalData
from typing import Union, List, Optional


def _inject_genomic_info_to_data(
    data: Union[MultimodalData, UnimodalData],
    genome: str,
) -> None:
    genomic_df = pd.read_csv(predefined_gene_orders[genome], sep="\t", header=None)
    genomic_df.columns = ["featureid", "chromosome", "start", "end"]
    genomic_df.set_index("featureid", inplace=True)

    var_df = data.var.reset_index().set_index("featureid")
    var_df = var_df.join(genomic_df)
    n_miss = genomic_df.shape[0] - var_df.shape[0]
    if n_miss > 0:
        logger.warning(f"InferCNV ignores {n_miss} genes which do not exist in data!")

    data.var = var_df.reset_index().set_index("featurekey")


def calc_infercnv(
    data: Union[MultimodalData, UnimodalData],
    genome: str,
    reference_key: str,
    reference_cat: List[str],
    mat_key: Optional[str] = None,
    window_size: int = 100,
    cnv_res: str = "cnv",
    chunk_size: int = 5000,
    n_jobs: int = -1,
) -> None:
    _inject_genomic_info_to_data(data, genome)

    # 1. Initial CNV scores
    # Step 1a. logTPM
    log_norm(data, norm_count=1e6, base_matrix=process_mat_key(data, mat_key), target_matrix="log_tpm")
    X = data.get_matrix("log_tpm")
    if issparse(X) and (not isinstance(X, csr_matrix)):
        X = X.tocsr()

    ## Step 1b. Center by gene means

    genomic_df = data.var.loc[:, ["chromosome", "start", "end"]]
    chromosomes = natsorted([x for x in genomic_df["chromosome"].unique() if x.startswith("chr")])
    chr2gene_idx = {}
    for chr in chromosomes:
        genes = genomic_df.loc[genomic_df["chromosome"]==chr].sort_values("start").index.values
        chr2gene_idx[chr] = genomic_df.index.get_indexer(genes)
    n_jobs = eff_n_jobs(n_jobs)
    chunks = process_map(
            _infercnv_chunk,
            [X[i:(i+chunk_size), :] for i in range(0, data.shape[0], chunk_size)],
            itertools.repeat(chr2gene_idx),
            itertools.repeat(window_size),
            tqdm_class=tqdm,
            max_workers=n_jobs,
        )
    res = scipy.sparse.vstack(chunks)

    ## 2. Final CNV scores

    data.obsm[f"X_{cnv_res}"] = res
    data.uns[cnv_res] = dict(zip(chromosomes, np.cumsum([0] + [idx.size for idx in chr2gene_idx.values()])))


def _infercnv_chunk(x, chr2gene_idx, window_size):
    from pegasus.cylib.fast_utils import calc_running_mean

    # Clip

    # Window smoothing
    running_means = []
    for chr in chr2gene_idx:
        tmp_x = x[:, chr2gene_idx[chr]].toarray()
        m, n = tmp_x.shape
        running_means.append(calc_running_mean(tmp_x, m, n, window_size))
    x_smoothed = np.hstack(running_means)

    # Center by cell medians

    x_res = x_smoothed
    return csr_matrix(x_res)
