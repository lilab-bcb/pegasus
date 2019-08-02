#!/usr/bin/env python

import time
import numpy as np
import pandas as pd
import os.path
from scipy.io import mmread
from scipy.sparse import csr_matrix
import tables

from typing import List
from . import Array2D, MemData

import anndata



def load_10x_h5_file_v2(h5_in: 'tables.File', fn: str, ngene: int = None) -> 'MemData':
	"""Load 10x v2 format matrix from hdf5 file
	
	Parameters
	----------

	h5_in : tables.File 
		An instance of tables.File class that is connected to a 10x v2 formatted hdf5 file.
	fn : `str`
		File name, can be used to indicate channel-specific name prefix.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.

	Returns
	-------
	
	An MemData object containing genome-Array2D pair per genome.

	Examples
	--------
	>>> io.load_10x_h5_file_v2(h5_in)
	"""	
	
	data = MemData()
	for group in h5_in.list_nodes("/", "Group"):
		genome = group._v_name

		M, N = h5_in.get_node("/" + genome + "/shape").read()
		mat = csr_matrix((h5_in.get_node("/" + genome + "/data").read(), h5_in.get_node("/" + genome + "/indices").read(), h5_in.get_node("/" + genome + "/indptr").read()), shape = (N, M))
		
		barcodes = h5_in.get_node("/" + genome + "/barcodes").read().astype(str)
		ids = h5_in.get_node("/" + genome + "/genes").read().astype(str)
		names = h5_in.get_node("/" + genome + "/gene_names").read().astype(str)

		array2d = Array2D({'barcodekey' : barcodes}, {'featurekey' : ids, 'featurename' : names}, mat)
		array2d.filter(ngene)
		array2d.separate_channels(fn)

		data.addData(genome, array2d)

	return data



def load_10x_h5_file_v3(h5_in: 'tables.File', fn: str, ngene: int = None) -> 'MemData':
	"""Load 10x v3 format matrix from hdf5 file
	
	Parameters
	----------

	h5_in : tables.File 
		An instance of tables.File class that is connected to a 10x v3 formatted hdf5 file.
	fn : `str`
		File name, can be used to indicate channel-specific name prefix.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.

	Returns
	-------
	
	An MemData object containing genome-Array2D pair per genome.

	Examples
	--------
	>>> io.load_10x_h5_file_v3(h5_in)
	"""	
	
	M, N = h5_in.get_node("/matrix/shape").read()
	bigmat = csr_matrix((h5_in.get_node("/matrix/data").read(), h5_in.get_node("/matrix/indices").read(), h5_in.get_node("/matrix/indptr").read()), shape = (N, M))
	barcodes = h5_in.get_node("/matrix/barcodes").read().astype(str)
	genomes = h5_in.get_node("/matrix/features/genome").read().astype(str)
	ids = h5_in.get_node("/matrix/features/id").read().astype(str)
	names = h5_in.get_node("/matrix/features/name").read().astype(str)

	data = MemData()
	for genome in np.unique(genomes):
		idx = genomes == genome

		barcode_metadata = {'barcodekey' : barcodes}
		feature_metadata = {'featurekey' : ids[idx], 'featurename' : names[idx]}
		mat = bigmat[:, idx].copy()
		array2d = Array2D(barcode_metadata, feature_metadata, mat)
		array2d.filter(ngene)
		array2d.separate_channels(fn)

		data.addData(genome, array2d)

	return data



def load_10x_h5_file(input_h5: str, ngene: int = None) -> 'MemData':
	"""Load 10x format matrix (either v2 or v3) from hdf5 file
	
	Parameters
	----------

	input_h5 : `str`
		The matrix in 10x v2 or v3 hdf5 format.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.

	Returns
	-------
	
	An MemData object containing genome-Array2D pair per genome.

	Examples
	--------
	>>> io.load_10x_h5_file('example_10x.h5')
	"""	
	
	fn = os.path.basename(input_h5)[:-3]

	data = None
	with tables.open_file(input_h5) as h5_in:
		try:
			node = h5_in.get_node("/matrix")
			data = load_10x_h5_file_v3(h5_in, fn, ngene)
		except tables.exceptions.NoSuchNodeError:
			data = load_10x_h5_file_v2(h5_in, fn, ngene)

	return data



def determine_file_name(path: str, names: List[str], errmsg: str) -> str:
	""" Try several file name options and determine which one is correct.
	"""
	for name in names:
		file_name = os.path.join(path, name)
		if os.path.exists(file_name):
			return file_name
	raise ValueError(errmsg)


def load_one_mtx_file(path: str) -> 'Array2D':
	"""Load one gene-count matrix in mtx format into an Array2D object
	"""
	mtx_file = determine_file_name(path, ['matrix.mtx.gz', 'matrix.mtx'], 'Expression matrix in mtx format is not found')
	mat = csr_matrix(mmread(mtx_file).T)

	barcode_file = determine_file_name(path, ['cells.tsv.gz', 'barcodes.tsv.gz', 'barcodes.tsv'], 'Barcode metadata information is not found')
	feature_file = determine_file_name(path, ['genes.tsv.gz', 'features.tsv.gz', 'genes.tsv'], 'Feature metadata information is not found')
	format_type = 'HCA DCP' if os.path.basename(barcode_file).startswith('cells') else ('10x v3' if feature_file.startswith('features') else '10x v2')

	if format_type == 'HCA DCP':
		barcode_metadata = pd.read_csv(barcode_file, sep = '\t', header = 0)
		assert 'cellkey' in barcode_metadata
		barcode_metadata.rename(columns = {'cellkey': 'barcodekey'}, inplace = True)

		feature_metadata = pd.read_csv(feature_file, sep = '\t', header = 0)
	elif format_type == '10x v3':
		barcode_metadata = pd.read_csv(barcode_file, sep = '\t', header = None, names = ['barcodekey'])
		feature_metadata = pd.read_csv(feature_file, sep = '\t', header = None, names = ['featurekey', 'featurename', 'featuretype'])
	else:
		barcode_metadata = pd.read_csv(barcode_file, sep = '\t', header = None, names = ['barcodekey'])
		feature_metadata = pd.read_csv(feature_file, sep = '\t', header = None, names = ['featurekey', 'featurename'])

	return Array2D(barcode_metadata, feature_metadata, mat)


def load_mtx_file(path: str, genome: str = None) -> 'MemData':
	"""Load gene-count matrix from Market Matrix files (10x v2, v3 and HCA DCP formats)
	
	Parameters
	----------

	path : `str`
		Path to mtx files. The directory impiled by path should either contain matrix, feature and barcode information, or folders containg these information.
	genome : `str`, optional (default: None)
		Genome name of the matrix. If None, genome will be inferred from path

	Returns
	-------
	
	An MemData object containing a genome-Array2D pair.

	Examples
	--------
	>>> io.load_10x_h5_file('example_10x.h5')
	"""	

	path = os.path.expanduser(os.path.expandvars(path))
	if not os.path.isdir(path):
		path = os.path.dirname(path)

	data = MemData()
	if os.path.exists(os.path.join(path, 'matrix.mtx.gz')) or os.path.exists(os.path.join(path, 'matrix.mtx')):
		if genome is None:
			genome = os.path.basename(path)
		data.addData(genome, load_one_mtx_file(path))
	else:
		for dir_entry in os.scandir(path):
			if dir_entry.is_dir():
				data.addData(dir_entry.name, load_one_mtx_file(dir_entry.path))

	return data



def load_dge_file(input_file: str, genome: str) -> 'MemData':
	"""Load DGE-format matrix from *dge.txt.gz, this is the format for Drop-seq, Seqwell and Celsee.

	Parameters
	----------

	input_file : `str`
		The matrix in dropseq format.
	genome : `str`
		The genome reference.

	Returns
	-------

	An MemData object containing a genome-Array2D pair.

	Examples
	--------
	>>> io.load_dge_file('example.umi.dge.txt.gz', 'GRCh38')
	"""

	df = pd.read_csv(input_file, header = 0, index_col = 0, sep = '\t', compression = 'gzip')
	mat = csr_matrix(df.values.T)
	barcode_metadata = {'barcodekey': df.columns.values}
	genes = df.index.values
	feature_metadata = {'featurekey': genes, 'featurename': genes}

	data = MemData()
	data.addData(genome, Array2D(barcode_metadata, feature_metadata, mat))

	return data



def load_csv_file(input_csv: str, genome: str) -> 'MemData':
	"""Load count matrix from a CSV file, including loading ADT matrix

	Parameters
	----------

	input_csv : `str`
		The CSV file containing the count matrix.
	genome : `str`
		The genome reference.

	Returns
	-------

	An MemData object containing a genome-Array2D pair.
	
	Examples
	--------
	>>> io.load_csv_file('example_ADT.csv')
	"""

	input_csv = os.path.expanduser(os.path.expandvars(input_csv))
	path = os.path.dirname(input_csv)
	base = os.path.basename(input_csv)

	converter = float if base.startswith('expression') else int

	barcodes = []
	names = []
	stacks = []
	with open(input_csv) as fin:
		barcodes = next(fin).strip().split(',')[1:]
		for line in fin:
			fields = line.strip().split(',')
			names.append(fields[0])
			stacks.append([converter(x) for x in fields[1:]])

	mat = csr_matrix(np.stack(stacks, axis = 1))
	barcode_metadata = {'barcodekey' : barcodes}
	feature_metadata = {'featurekey' : names, 'featurename' : names}

	if base == 'expression.csv':
		barcode_file = os.path.join(path, 'cells.csv')
		if os.path.exists(barcode_file):
			barcode_metadata = pd.read_csv(barcode_file, sep = ',', header = 0)
			assert 'cellkey' in barcode_metadata
			barcode_metadata.rename(columns = {'cellkey': 'barcodekey'}, inplace = True)
		
		feature_file = os.path.join(path, 'genes.csv')
		if os.path.exists(feature_file):
			feature_metadata = pd.read_csv(feature_file, sep = ',', header = 0)

	data = MemData()
	data.addData(genome, Array2D(barcode_metadata, feature_metadata, mat))

	return data



def load_scCloud_h5_file(input_h5: str, ngene: int = None, select_singlets: bool = False) -> 'MemData':
	"""Load matrices from scCloud-format hdf5 file
	
	Parameters
	----------

	input_h5 : `str`
		scCloud-format hdf5 file.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.
	select_singlets: `bool`, optional (default: False)
		If only load singlets.

	Returns
	-------
	
	An MemData object containing genome-Array2D pair per genome.

	Examples
	--------
	>>> io.load_scCloud_h5_file('example.scCloud.h5')
	"""	

	cite_seq_name = None
	selected_barcodes = None

	data = MemData()
	with tables.open_file(input_h5) as h5_in:
		for group in h5_in.list_nodes("/", "Group"):
			genome = group._v_name

			M, N = h5_in.get_node("/" + genome + "/shape").read()
			mat = csr_matrix((h5_in.get_node("/" + genome + "/data").read(), h5_in.get_node("/" + genome + "/indices").read(), \
							h5_in.get_node("/" + genome + "/indptr").read()), shape = (N, M))

			barcode_metadata = {}
			for node in h5_in.walk_nodes("/" + genome + "/_barcodes", "Array"):
				values = node.read()
				if values.dtype.kind == 'S':
					values = values.astype(str)
				barcode_metadata[node.name] = values

			feature_metadata = {}
			for node in h5_in.walk_nodes("/" + genome + "/_features", "Array"):
				values = node.read()
				if values.dtype.kind == 'S':
					values = values.astype(str)
				feature_metadata[node.name] = values

			array2d = Array2D(barcode_metadata, feature_metadata, mat)
			if genome.startswith("CITE_Seq"):
				cite_seq_name = genome
			else:
				array2d.filter(ngene, select_singlets)
				selected_barcodes = array2d.get_metadata('barcodekey')

			data.addData(genome, array2d)

	if (cite_seq_name is not None) and (selected_barcodes is not None):
		array2d = data.getData(cite_seq_name)
		selected = array2d.get_metadata('barcodekey').isin(selected_barcodes)
		array2d.trim(selected)

	return data



def read_input(input_file: str, return_type = 'AnnData', genome: str = None, concat_matrices: bool = False, h5ad_mode: str = 'r+', select_singlets: bool = False) -> 'MemData or AnnData or List[AnnData]':
	"""Load data into memory.

	This function is used to load input data into memory. Inputs can be in 10x genomics v2 & v3 formats (hdf5 or mtx), HCA DCP mtx and csv formats, Drop-seq dge format, and CSV format.

	Parameters
	----------

	input_file : `str`
		Input file name.
	return_type : `str`
		Return object type, can be either 'MemData' or 'AnnData'.
	genome : `str`, optional (default: None)
		A string contains comma-separated genome names. scCloud will read all matrices matching the genome names. If genomes is None, all matrices will be considered.
	concat_matrices : `boolean`, optional (default: False)
		If input file contains multiple matrices, if concatenate them into one AnnData object or return a list of AnnData objects.
	h5ad_mode : `str`, optional (default: `r+`)
		If input is in h5ad format, the backed mode for loading the data. mode could be 'a', 'r', 'r+'. 'a' refers to load all into memory.
	select_singlets : `bool`, optional (default: False)
		If only keep DemuxEM-predicted singlets when loading data.

	Returns
	-------
	`MemData` object or `anndata` object or a list of `anndata` objects
		An `MemData` object or `anndata` object or a list of `anndata` objects containing the count matrices.

	Examples
	--------
	>>> adata = io.read_input('example_10x.h5', genomes = 'mm10')
	>>> adata = io.read_input('example.h5ad', mode = 'r+')
	>>> adata = io.read_input('example_ADT.csv')
	"""

	start = time.time()

	if input_file.endswith('.scCloud.h5'):
		data = load_scCloud_h5_file(input_file, select_singlets = select_singlets)
	elif input_file.endswith('.h5'):
		data = load_10x_h5_file(input_file)
	elif input_file.endswith('.h5ad'):
		data = anndata.read_h5ad(input_file, backed = (False if h5ad_mode == 'a' else h5ad_mode))
	elif input_file.endswith('.mtx') or input_file.endswith('.mtx.gz') or os.path.isdir(input_file):
		data = load_mtx_file(input_file, genome)
	elif input_file.endswith('dge.txt.gz'):
		assert genome is not None
		data = load_dge_file(input_file, genome)
	elif input_file.endswith('.csv'):
		assert genome is not None
		data = load_csv_file(input_file, genome)
	else:
		raise ValueError("Unrecognized file type")

	data.restrain_keywords(genome)

	if return_type == 'AnnData':
		data = data.convert_to_anndata(concat_matrices = concat_matrices)

	end = time.time()
	print("Read input is finished. Time spent = {:.2f}s.".format(end - start))

	return data



def write_output(data: 'MemData or AnnData', output_name: str) -> None:
	""" Write data back to disk.

	This function is used to write data back to disk. 

	Parameters
	----------

	data : `MemData` or `AnnData`
		data to write back, can be either an MemData or AnnData object.
	output_name : `str`
		output file name. MemData ends with suffix '.scCloud.h5' and AnnData ends with suffix '.h5ad'

	Returns
	-------
	None

	Examples
	--------
	>>> io.write_output(adata, 'test')
	"""

	start = time.time()

	if isisntance(data, MemData):
		data.write_h5_file(output_name + '.scCloud.h5')
	else:
		data.write(output_name + '.h5ad')

	end = time.time()
	print("Write output is finished. Time spent = {:.2f}s.".format(end - start))
