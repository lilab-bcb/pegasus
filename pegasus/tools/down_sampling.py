import time
import numpy as np
import pandas as pd
import tables
import numba
from numba import njit
from scipy.sparse import csr_matrix

from pegasus.io import write_output, Array2D, MemData


def sample_hypergeom(count, total_reads, n_sample, random_state):
    np.random.seed(random_state)
    new_count = np.zeros(count.size, dtype=int)
    for i in range(count.size):
        new_count[i] = np.random.hypergeometric(
            count[i], total_reads - count[i], n_sample
        )
        total_reads -= count[i]
        n_sample -= new_count[i]
    return new_count


@njit
def generate_sparse_matrix(barcode_idx, feature_idx):
    row_ind = []
    col_ind = []
    data = []

    start = 0
    for i in range(1, barcode_idx.size):
        if barcode_idx[start] != barcode_idx[i] or feature_idx[start] != feature_idx[i]:
            row_ind.append(barcode_idx[start])
            col_ind.append(feature_idx[start])
            data.append(i - start)
            start = i

    row_ind.append(barcode_idx[start])
    col_ind.append(feature_idx[start])
    data.append(barcode_idx.size - start)

    return row_ind, col_ind, data


def down_sample(molecule_info_file, output_name, total_reads, n_sample, random_state=0):
    with tables.open_file(molecule_info_file) as h5_in:
        barcode_idx = h5_in.get_node("/barcode_idx").read()
        feature_idx = h5_in.get_node("/feature_idx").read()
        count = h5_in.get_node("/count").read()

        new_count = sample_hypergeom(count, total_reads, n_sample, random_state)

        idx = new_count > 0
        barcode_idx = barcode_idx[idx]
        feature_idx = feature_idx[idx]

        barcodes = h5_in.get_node("/barcodes").read()
        gene_ids = h5_in.get_node("/features/id").read()
        gene_names = h5_in.get_node("/features/name").read()

        genome = h5_in.get_node("/barcode_info/genomes").read()[0].decode()

        row_ind, col_ind, data = generate_sparse_matrix(barcode_idx, feature_idx)
        mat = csr_matrix(
            (data, (row_ind, col_ind)), shape=(barcodes.size, gene_ids.size)
        )

        data = MemData()
        data.addData(
            genome,
            Array2D(
                {"barcodekey": barcodes},
                {"featurekey": gene_ids, "featurename": gene_names},
                mat,
            ),
        )

        write_output(data, output_name)

        print("Subsampled raw matrix is generated!")
