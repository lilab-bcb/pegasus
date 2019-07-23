import numpy as np
import pandas as pd
import time
from scipy.sparse import issparse
from . import read_input

def scp_write_coords(data, output_name):
	cluster_labels = []
	for col_name in data.obs.columns:
		if col_name.find('labels') >= 0:
			cluster_labels.append(col_name)
	df_labels = data.obs[cluster_labels]
	clu_str = group_str = ''
	if len(cluster_labels) > 0:
		clu_str = ''.join(['\t' + x for x in cluster_labels])
		group_str = ''.join(['\tgroup'] * len(cluster_labels))
	basis_set = set(data.obsm.dtype.names)
	for basis in ['X_tsne', 'X_fitsne', 'X_umap', 'X_diffmap_pca', 'X_fle']:
		if basis in basis_set:
			coords = ['X', 'Y'] if basis != 'X_diffmap_pca' else ['X', 'Y', 'Z']
			coo_str = '\t'.join(coords)
			num_str = '\t'.join(['numeric'] * len(coords))
			coord_file = "{}.scp.{}.coords.txt".format(output_name, basis)
			with open(coord_file, 'w') as fout:
				fout.write("NAME\t{coo}{clu}\n".format(coo = coo_str, clu = clu_str))
				fout.write("TYPE\t{coo}{clu}\n".format(coo = num_str, clu = group_str))
			df_out = pd.DataFrame(data.obsm[basis][:,0:len(coords)], columns = coords, index = data.obs_names)
			df_out = pd.concat([df_out, df_labels], axis = 1)
			df_out.to_csv(coord_file, sep = "\t", header = False, mode = 'a')
			print("Coordinate file {} is written.".format(coord_file))

def scp_write_metadata(data, output_name):
	ban = ['n_genes', 'n_counts', 'percent_mito', 'pseudotime']
	meta = []
	for col_name in data.obs.columns:
		if (col_name not in ban) and (col_name.find('labels') < 0):
			meta.append(col_name)
	meta_str = ''.join(['\t' + x for x in meta])
	group_str = ''.join(['\tgroup'] * len(meta))
	metadata_file = "{}.scp.metadata.txt".format(output_name)
	with open(metadata_file, 'w') as fout:
		fout.write("NAME{meta}\n".format(meta = meta_str))
		fout.write("TYPE{meta}\n".format(meta = group_str))
	data.obs[meta].to_csv(metadata_file, sep = '\t', header = False, mode = 'a')
	print("Metadata file {} is written.".format(metadata_file))


def write_market_matrix(mtx_file, X, round_to):
	with open(mtx_file, 'w') as fout:
		fmt_str = "{{}} {{}} {{:.{}f}}\n".format(round_to)
		fout.write("%%MatrixMarket matrix coordinate real general\n%\n")
		if issparse(X):
			X = X.tocoo()
			fout.write("{} {} {}\n".format(X.shape[0], X.shape[1], (X.data != 0).sum()))
			for x, y, value in zip(X.row, X.col, X.data):
				if value != 0:
					fout.write(fmt_str.format(x + 1, y + 1, value))
		else:
			fout.write("{} {} {}\n".format(X.shape[0], X.shape[1], np.count_nonzero(X)))
			for x, y in zip(*X.nonzero()):
				fout.write(fmt_str.format(x + 1, y + 1, X[x, y]))

def write_dense_matrix(expr_file, data, round_to):
	fmt_str = "{{:.{}f}}".format(round_to)
	exprmat = data.X.transpose()
	if issparse(exprmat):
		exprmat = exprmat.toarray()
	with open(expr_file, 'w') as fout:
		fout.write("GENE\t" + "\t".join(data.obs_names) + "\n")
		for i in range(exprmat.shape[0]):
			fout.write(data.var_names[i] + "\t" + "\t".join([fmt_str.format(x) if x > 0.0 else "0" for x in exprmat[i,:]]) + "\n")
	print("Expression matrix {} is written.".format(expr_file))


def scp_write_expression(data, output_name, is_sparse = True, round_to = 2):
	if is_sparse:
		barcode_file = "{}.scp.barcodes.tsv".format(output_name)
		with open(barcode_file, 'w') as fout:
			fout.write('\n'.join(data.obs_names) + '\n')
		print("Barcode file {} is written.".format(barcode_file))
		gene_file = "{}.scp.genes.tsv".format(output_name)
		df = pd.DataFrame({'gene_names' : data.var_names, 'gene_ids' : data.var['gene_ids']})[['gene_ids', 'gene_names']]
		with open(gene_file, 'w') as fout:
			df.to_csv(fout, sep = ' ', header = False, index = False)
		print("Gene file {} is written.".format(gene_file))
		mtx_file = "{}.scp.matrix.mtx".format(output_name)
		write_market_matrix(mtx_file, data.X.transpose().tocsr(), round_to)
		print("Matrix file {} is written.".format(mtx_file))
	else:
		expr_file = "{}.scp.expr.txt".format(output_name)
		write_dense_matrix(expr_file, data, round_to)

def run_scp_output(input_h5ad_file, output_name, is_sparse, round_to):
	adata = read_input(input_h5ad_file, mode = 'a')
	start = time.time()
	scp_write_coords(adata, output_name)
	scp_write_metadata(adata, output_name)
	scp_write_expression(adata, output_name, is_sparse, round_to)
	end = time.time()
	print("Time spent for generating SCP outputs is {:.2f}s.".format(end - start))
