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
from pegasusio import MultimodalData, UnimodalData, timer
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
        logger.warning(f"InferCNV ignores {n_miss} genes which do not exist in data.")

    data.var = var_df.reset_index().set_index("featurekey")


@timer(logger=logger)
def calc_infercnv(
    data: Union[MultimodalData, UnimodalData],
    genome: str,
    reference_key: Optional[str] = None,
    reference_cat: Optional[List[str]] = None,
    mat_key: Optional[str] = None,
    run_log_norm: bool = True,
    noise: float = 0.2,
    window_size: int = 100,
    res_key: str = "cnv",
    exclude_chromosomes: Optional[List[str]] = ["chrM", "chrX", "chrY"],
    chunk_size: int = 5000,
    n_jobs: int = -1,
) -> None:
    var_genomic_cols = [col for col in data.var.columns if col in ["chromosome", "start", "end"]]
    if len(var_genomic_cols) < 3:
        _inject_genomic_info_to_data(data, genome)

    # I. Initial CNV scores
    # log-norm TP100K
    source_mat_key = mat_key
    if run_log_norm:
        mat_key = process_mat_key(data, mat_key)
        source_mat_key = f"{mat_key}.log_norm"
        log_norm(data, base_matrix=mat_key, target_matrix=source_mat_key)

    # Exclude genes with no chromosome position or on chromosomes not included
    var_mask = data.var["chromosome"].isnull()
    n_no_pos = np.sum(var_mask)
    if n_no_pos > 0:
        logger.warning(f"Skip {n_no_pos} genes which have no genomic position annotated.")
    if exclude_chromosomes:
        logger.warning(f"Exclude chromosomes: {exclude_chromosomes}.")
        var_mask = var_mask | data.var["chromosome"].isin(exclude_chromosomes)
    data_sub = data[:, ~var_mask].copy()

    X = data_sub.get_matrix(source_mat_key)
    if issparse(X) and (not isinstance(X, csr_matrix)):
        X = X.tocsr()

    genomic_df = data_sub.var.loc[:, ["chromosome", "start", "end"]]
    chromosomes = natsorted([x for x in genomic_df["chromosome"].unique() if x.startswith("chr")])
    chr2gene_idx = {}
    for chr in chromosomes:
        genes = genomic_df.loc[genomic_df["chromosome"]==chr].sort_values("start").index.values
        chr2gene_idx[chr] = genomic_df.index.get_indexer(genes)

    n_jobs = eff_n_jobs(n_jobs)

    logger.info("Calculate initial CNV scores")
    chunks = process_map(
            _infercnv_chunk,
            [X[i:(i+chunk_size), :] for i in range(0, data.shape[0], chunk_size)],
            itertools.repeat(chr2gene_idx),
            itertools.repeat(window_size),
            tqdm_class=tqdm,
            max_workers=n_jobs,
        )
    cnv_init = np.vstack(chunks)

    ## II. Final CNV scores
    if reference_cat and reference_key:
        base_means = _get_reference_means(data_sub, cnv_init, reference_key, reference_cat)
        base_min = np.min(base_means, axis=0)
        base_max = np.max(base_means, axis=0)
        #cnv_final = np.apply_along_axis(lambda row: _correct_by_reference_per_cell(row, base_min, base_max, noise), axis=1, arr=cnv_init)
        logger.info("Calculate final CNV scores")
        chunks = process_map(
            _correct_by_reference_per_cell,
            [cnv_init[i:(i+chunk_size), :] for i in range(0, data.shape[0], chunk_size)],
            itertools.repeat(base_min),
            itertools.repeat(base_max),
            itertools.repeat(noise),
            tqdm_class=tqdm,
            max_workers=n_jobs,
        )
        cnv_final = np.vstack(chunks)
    else:
        logger.info("No reference category is specified. The initial CNV scores are already the final scores.")
        cnv_final = cnv_init

    data.obsm[f"X_{res_key}"] = csr_matrix(cnv_final)
    data.uns[res_key] = dict(zip(chromosomes, np.cumsum([0] + [idx.size for idx in chr2gene_idx.values()])))


def _infercnv_chunk(x, chr2gene_idx, window_size):
    from pegasus.cylib.fast_utils import calc_running_mean

    # Step 1. Window smoothing
    running_means = []
    for chr in chr2gene_idx:
        tmp_x = x[:, chr2gene_idx[chr]].toarray()
        m, n = tmp_x.shape
        if n < window_size:
            window_size = n
        running_means.append(calc_running_mean(tmp_x, m, n, window_size))
    x_smoothed = np.hstack(running_means)

    # Step 2. Center by cell medians
    x_cell_centered = x_smoothed - np.median(x_smoothed, axis=1)[:, np.newaxis]

    return x_cell_centered


def _get_reference_means(data, cnv_init, reference_key, reference_cat):
    if isinstance(reference_cat, str):
        reference_cat = [reference_cat]
    reference_cat = np.array(reference_cat)
    is_observed = np.isin(reference_cat, data.obs[reference_key])
    if not np.all(is_observed):
        ref_missing = reference_cat[~is_observed]
        logger.warning(f"Omit reference categories missing in data: {ref_missing}!")
        reference_cat = reference_cat[is_observed]

    base_mean_list = []
    for cat in reference_cat:
        idx_cat = data.obs.index.get_indexer(data.obs.loc[data.obs[reference_key]==cat].index.values)
        base_mean_list.append(np.mean(cnv_init[idx_cat, :], axis=0))
    return np.vstack(base_mean_list)


def _correct_by_reference_per_cell(cnv_init_row, base_min, base_max, noise):
    cnv_final_row = np.zeros(cnv_init_row.shape, dtype=cnv_init_row.dtype)
    above_max = cnv_init_row > (base_max + noise)
    cnv_final_row[above_max] = (cnv_init_row - base_max)[above_max]
    below_min = cnv_init_row < (base_min - noise)
    cnv_final_row[below_min] = (cnv_init_row - base_min)[below_min]
    return cnv_final_row
