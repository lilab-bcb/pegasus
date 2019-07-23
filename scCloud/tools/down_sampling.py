import time
import numpy as np 
import pandas as pd
import tables
import numba
from numba import njit
from scipy.sparse import csr_matrix

from . import write_10x_h5_file


def sample_hypergeom(count, total_reads, n_sample, random_state):
	np.random.seed(random_state)
	new_count = np.zeros(count.size, dtype = int)
	for i in range(count.size):
		new_count[i] = np.random.hypergeometric(count[i], total_reads - count[i], n_sample)
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

def down_sample(molecule_info_file, output_raw_file, total_reads, n_sample, random_state = 0):
	with tables.open_file(molecule_info_file) as h5_in:
		barcode_idx = h5_in.get_node('/barcode_idx').read()
		feature_idx = h5_in.get_node('/feature_idx').read()
		count = h5_in.get_node('/count').read()

		new_count = sample_hypergeom(count, total_reads, n_sample, random_state)

		idx = new_count > 0
		barcode_idx = barcode_idx[idx]
		feature_idx = feature_idx[idx]

		barcodes = h5_in.get_node('/barcodes').read()
		gene_ids = h5_in.get_node('/features/id').read()
		gene_names = h5_in.get_node('/features/name').read()

		channel = {}

		row_ind, col_ind, data = generate_sparse_matrix(barcode_idx, feature_idx)
		channel["matrix"] = csr_matrix((data, (row_ind, col_ind)), shape = (barcodes.size, gene_ids.size))
		channel["barcodes"] = barcodes
		channel["genes"] = gene_ids
		channel["gene_names"] = gene_names

		genome = h5_in.get_node('/barcode_info/genomes').read()[0].decode()

		results = {genome : channel}
		write_10x_h5_file(output_raw_file, results)
		print("Subsampled raw matrix is generated!")

