import time
import os.path
import numpy as np
import anndata
from scipy.sparse import csr_matrix, hstack
from . import col_attrs, excluded, load_10x_h5_file, load_dropseq_file
import pandas as pd
import scipy


def read_10x_h5_file(input_h5, genome = None, return_a_dict = False, select_singlets = False):
	"""Load 10x-format matrices from the h5 file into a series of h5ad objects

	Parameters
	----------

	input_h5 : `str`
		The matricies in h5 format.
	genome : `str`, optional (default: None)
		A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
	return_a_dict : `boolean`, optional (default: False)
		If input file contains multiple genome groups, if concatenate them into one h5ad object or return a dictionary of genome-h5ad pairs. If this option is on, return a dict.
	select_singlets : `int`, optional (default: False)
		If only keep DemuxEM-predicted singlets when loading data.

	Returns
	-------

	`anndata` object or a dictionary of `anndata` objects
		An `anndata` object or a dictionary of `anndata` objects containing the count matrices.

	Examples
	--------
	>>> tools.read_10x_h5_file('example_10x.h5')
	"""

	gdmap = load_10x_h5_file(input_h5, select_singlets = select_singlets) # gdmap , genome-data map
	if genome is not None: # remove genomes not in the list
		requested_genomes = set(genome.split(','))
		available_genomes = set(gdmap)
		remove_set = available_genomes - requested_genomes
		for gname in remove_set:
			gdmap.pop(gname)
		for genome in requested_genomes:
			if genome not in gdmap:
				raise ValueError('Genome {} does not exist in {}.'.format(genome, ', '.join(available_genomes)))

	results = {}
	for genome, data in gdmap.items():
		# obs_dict
		barcodes = data["barcodes"].astype(str)
		if np.chararray.endswith(barcodes, '-1').sum() == barcodes.size:
			barcodes = np.vectorize(lambda x: x[:-2])(barcodes)

		def extract_channel(barcode):
			fields = barcode.split('-')
			if len(fields) == 2 and fields[-1].isdigit():
				return fields[-1]
			nshift = 2 if fields[-1] == '1' else 1
			return '-'.join(fields[:-nshift])

		obs_dict = {"obs_names" : barcodes, "Channel" : np.vectorize(extract_channel)(barcodes)}

		for attr, value in data.items():
			if (attr not in col_attrs) and (attr not in excluded):
				obs_dict[attr] = value.astype(str)
		# var_dict
		var_dict = {"var_names" : (data["gene_names"] if "gene_names" in data else data["antibody_names"]).astype(str)}
		if "genes" in data:
			var_dict["gene_ids"] = data["genes"].astype(str)

		# convert matrix to float32 if necessary
		mat = data["matrix"]
		if mat.dtype == np.int32:
			mat.dtype = np.float32
			orig_dat = mat.data.view(np.int32)
			mat.data[:] = orig_dat

		# construct h5ad object
		results[genome] = anndata.AnnData(X = mat, obs = obs_dict, var = var_dict)
		results[genome].uns["genome"] = genome

	if len(results) == 1:
		results = results[next(iter(results))]
	elif not return_a_dict:
		Xs = [] # a list of csr matrices
		vn_vec = [] # var_names vec
		gi_vec = [] # gene_ids vec
		genomes = []
		for data in results.values():
			Xs.append(data.X)
			vn_vec.append(data.var_names.values)
			if "gene_ids" in data.var:
				gi_vec.append(data.var["gene_ids"].values)
			genomes.append(data.uns["genome"])
		var_dict = {"var_names" : np.concatenate(vn_vec)}
		if len(gi_vec) > 0:
			var_dict["gene_ids"] = np.concatenate(gi_vec)
		results = anndata.AnnData(X = hstack(Xs, format = 'csr'), obs = obs_dict, var = var_dict)
		results.uns["genome"] = ",".join(genomes)

	return results



def read_dropseq_file(input_file, genome):
	"""Load dropseq-format matrix from the dropseq file into a h5ad object

	Parameters
	----------

	input_file : `str`
		The matrix in dropseq format.
	genome : `str`
		The genome reference.

	Returns
	-------

	`anndata` object
		An `anndata` object containing the gene-count matrix

	Examples
	--------
	>>> tools.read_dropseq_file('example.umi.dge.txt.gz', 'GRCh38')
	"""

	results = load_dropseq_file(input_file, genome)
	df = results[genome]["matrix"]
	data = anndata.AnnData(X = csr_matrix(df.values.transpose()), obs = {"obs_names" : df.columns.values}, var = {"var_names" : df.index.values, "genes" : df.index.values})
	data.obs['Channel'] = input_file[:-len(".dge.txt.gz")]
	data.uns['genome'] = genome

	return data



def read_antibody_csv(input_csv):
	"""Load an ADT matrix from the csv file

	Parameters
	----------

	input_csv : `str`
		The CSV file containing ADT counts.

	Returns
	-------

	`anndata` object
		An `anndata` object containing the ADT count matrix.

	Examples
	--------
	>>> tools.read_antibody_csv('example_ADT.csv')
	"""

	barcodes = []
	antibody_names = []
	stacks = []
	with open(input_csv) as fin:
		barcodes = next(fin).strip().split(',')[1:]
		for line in fin:
			fields = line.strip().split(',')
			antibody_names.append(fields[0])
			stacks.append([int(x) for x in fields[1:]])
	data = anndata.AnnData(X = csr_matrix(np.stack(stacks, axis = 1)), obs = {"obs_names" : barcodes}, var = {"var_names" : antibody_names})

	return data



def read_input(input_file, genome = None, return_a_dict = False, mode = 'r+', select_singlets = False):
	"""Load data into memory.

	This function is used to load input data into memory. Inputs can be in .h5, .h5ad, .dge.txt.gz or .csv format.

	Parameters
	----------

	input_file : `str`
		Input file name.
	genome : `str`, .h5, optional (default: None)
		A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
	return_a_dict : `boolean`, .h5, optional (default: False)
		If input file contains multiple genome groups, if concatenate them into one h5ad object or return a dictionary of genome-h5ad pairs. If this option is on, return a dict.
	mode : `str`, .h5ad, optional (default: `r+`)
		If input is in h5ad format, the backed mode for loading the data. mode could be 'a', 'r', 'r+'. 'a' refers to load all into memory.
	select_singlets : `int`, optional (default: False)
		If only keep DemuxEM-predicted singlets when loading data.

	Returns
	-------
	`anndata` object or a dictionary of `anndata` objects
		An `anndata` object or a dictionary of `anndata` objects containing the count matrices.

	Examples
	--------
	>>> adata = tools.read_input('example_10x.h5', genome = 'mm10')
	>>> adata = tools.read_input('example.h5ad', mode = 'r+')
	>>> adata = tools.read_input('example_ADT.csv')
	"""

	start = time.time()

	if input_file.endswith('.h5'):
		data = read_10x_h5_file(input_file, genome, return_a_dict, select_singlets)
	elif input_file.endswith('.dge.txt.gz'):
		data = read_dropseq_file(input_file, genome)
	elif input_file.endswith('.h5ad'):
		data = anndata.read_h5ad(input_file, backed = (False if mode == 'a' else mode))
		# input_prefix = input_file[:-5]
		# index_file = input_prefix + '.knn.hnsw.bin'
		# if os.path.isfile(index_file):
		# 	import hnswlib
		# 	knn_index = hnswlib.Index(space = 'l2', dim = data.uns['knn_dim'])
		# 	knn_index.load_index(index_file)
		# 	data.uns['knn'] = knn_index

		# 	index_file = input_prefix + '.diffmap_knn.hnsw.bin'
		# 	assert os.path.isfile(index_file)
		# 	knn_index = hnswlib.Index(space = 'l2', dim = data.uns['diffmap_knn_dim'])
		# 	knn_index.load_index(index_file)
		# 	data.uns['diffmap_knn'] = knn_index
	elif input_file.endswith('.csv'):
		data = read_antibody_csv(input_file)
	elif input_file.endswith('.mtx') or input_file.endswith('.mtx.gz') or os.path.isdir(input_file):
		data = __read_mtx(input_file)
	else:
		raise ValueError("Unrecognized file type")


	end = time.time()
	print("Read input is finished. Time spent = {:.2f}s.".format(end - start))

	return data


def __read_mtx(path):
	if os.path.isdir(path):
		if os.path.exists(os.path.join(path, 'matrix.mtx')):
			path = os.path.join(path, 'matrix.mtx')
		else:
			path = os.path.join(path, 'matrix.mtx.gz')
		if not os.path.exists(path):
			raise ValueError('matrix.mtx not found')
	parent_directory = os.path.dirname(path)
	x = scipy.io.mmread(path)
	x = scipy.sparse.csr_matrix(x.T)

	is_10x_v3 = os.path.exists(os.path.join(parent_directory, 'features.tsv.gz'))
	# feature ID, name, type for v3 features.tsv.gz

	var = pd.read_csv(os.path.join(parent_directory, 'features.tsv.gz' if is_10x_v3 else 'genes.tsv'), sep='\t',
		header=None, index_col=1, names=['gene_ids', 'var_names', 'feature_type'] if is_10x_v3 else ['gene_ids', 'var_names'])
	var.index = anndata.utils.make_index_unique(var.index)
	obs = pd.read_csv(os.path.join(parent_directory, 'barcodes.tsv.gz' if is_10x_v3 else 'barcodes.tsv'), sep='\t',
		header=None, index_col=0, names=['obs_names'])
	cell_count, gene_count = x.shape
	if len(obs)!=cell_count:
		raise ValueError(
			"Wrong number of cells : matrix has {} cells, barcodes file has {}".format(cell_count, len(obs)))
	if len(var)!=gene_count:
		raise ValueError("Wrong number of genes : matrix has {} genes, genes file has {}".format(gene_count, len(var)))

	return anndata.AnnData(X=x, obs=obs, var=var)

def write_output(data, output_name):
	start = time.time()

	# if 'knn' in data.uns:
	# 	data.uns['knn'].save_index(output_name + '.knn.hnsw.bin')
	# 	del data.uns['knn']
	# if 'diffmap_knn' in data.uns:
	# 	data.uns['diffmap_knn'].save_index(output_name + '.diffmap_knn.hnsw.bin')
	# 	del data.uns['diffmap_knn']
	data.write(output_name + ".h5ad")

	end = time.time()
	print("Write main output is finished. Time spent = {:.2f}s.".format(end - start))
