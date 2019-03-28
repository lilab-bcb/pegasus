#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from scipy.sparse import csr_matrix, vstack
import tables
import copy
from subprocess import check_call

from . import col_attrs, excluded



def load_10x_h5_file_v2(h5_in, ngene = None, select_singlets = False):
	"""Load 10x v2 format matrix from hdf5 file
	
	Parameters
	----------

	h5_in : tables.File 
		An instance of tables.File class that is connected to a 10x v2 formatted hdf5 file.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.

	Returns
	-------
	
	A dictionary containing genome-channel pair per each genome. The channels are also dictionaries containing the count matricies in nsample x ngene format.

	Examples
	--------
	>>> tools.load_10x_h5_file_v2(h5_in)
	"""	
	
	def channel_select(channel, selected):
		channel["matrix"] =  channel["matrix"][selected, :]
		channel["barcodes"] = channel["barcodes"][selected]
		for attr in channel:
			if (attr not in col_attrs) and (attr not in excluded):
				channel[attr] = channel[attr][selected]
		return channel["barcodes"]



	results = {}

	cite_seq_name = None
	selected_barcodes = None

	for i, group in enumerate(h5_in.walk_groups()):
		if i > 0:
			genome = group._v_name
			
			channel = {}
			for node in h5_in.walk_nodes("/" + genome, "Array"):
				channel[node.name] = node.read()
			
			M, N = channel["shape"]
			channel["matrix"] = csr_matrix((channel["data"], channel["indices"], channel["indptr"]), shape = (N, M))
			channel.pop("data")
			channel.pop("indices")
			channel.pop("indptr")
			channel.pop("shape")


			if genome.startswith("CITE_Seq"):
				cite_seq_name = genome
			else:
				if (ngene is not None) or (select_singlets and "demux_type" in channel):
					selected = np.ones(N, dtype = bool)
					if ngene is not None:
						selected = selected & (channel["matrix"].getnnz(axis = 1) >= ngene)
					if select_singlets:
						selected = selected & (channel["demux_type"] == "singlet".encode())
						channel.pop("demux_type")
					selected_barcodes = channel_select(channel, selected)

			results[genome] = channel

	if (cite_seq_name is not None) and (selected_barcodes is not None):
		selected_barcodes = pd.Index(selected_barcodes)
		channel = results[cite_seq_name]
		citeseq_barcodes = pd.Index(channel["barcodes"])
		selected = citeseq_barcodes.isin(selected_barcodes)
		channel_select(channel, selected)

	return results



def load_10x_h5_file_v3(h5_in, ngene = None):
	"""Load 10x v3 format matrix from hdf5 file
	
	Parameters
	----------

	h5_in : tables.File 
		An instance of tables.File class that is connected to a 10x v3 formatted hdf5 file.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.

	Returns
	-------
	
	A dictionary containing genome-channel pair per each genome. The channels are also dictionaries containing the count matricies.

	Examples
	--------
	>>> tools.load_10x_h5_file_v3(h5_in)
	"""	
	
	M, N = h5_in.get_node("/matrix/shape").read()
	bigmat = csr_matrix((h5_in.get_node("/matrix/data").read(), h5_in.get_node("/matrix/indices").read(), h5_in.get_node("/matrix/indptr").read()), shape = (N, M))
	barcodes = h5_in.get_node("/matrix/barcodes").read()
	genomes = h5_in.get_node("/matrix/features/genome").read()
	ids = h5_in.get_node("/matrix/features/id").read()
	names = h5_in.get_node("/matrix/features/name").read()

	results = {}
	for genome in np.unique(genomes):
		idx = genomes == genome
		channel = {"genes" : ids[idx], "gene_names" : names[idx]}
		mat = bigmat[:, idx].copy()

		if ngene is not None:
			selected = mat.getnnz(axis = 1) >= ngene
			channel["matrix"] = mat[selected, :]
			channel["barcodes"] = barcodes[selected]
		else:
			channel["matrix"] = mat
			channel["barcodes"] = barcodes

		results[genome.decode()] = channel

	return results



def load_10x_h5_file(input_h5, ngene = None, select_singlets = False):
	"""Load 10x format matrix (either v2 or v3) from hdf5 file
	
	Parameters
	----------

	input_h5 : `str`
		The matrix in 10x v2 or v3 hdf5 format.
	ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode. Default is to keep all barcodes.

	Returns
	-------
	
	A dictionary containing genome-channel pair per each genome. The channels are also dictionaries containing the count matricies.

	Examples
	--------
	>>> tools.load_10x_h5_file('example_10x.h5')
	"""	
	
	results = None
	with tables.open_file(input_h5) as h5_in:
		try:
			node = h5_in.get_node("/matrix")
			results = load_10x_h5_file_v3(h5_in, ngene)
		except tables.exceptions.NoSuchNodeError:
			results = load_10x_h5_file_v2(h5_in, ngene, select_singlets)
	return results



def load_dropseq_file(input_file, genome):
	"""Load dropseq-format matrix from the dropseq file

	Parameters
	----------

	input_file : `str`
		The matrix in dropseq format.
	genome : `str`
		The genome reference.

	Returns
	-------

	A dictionary containing one genome-channel pair. The genome is the provided genome name. The channel is a dictionary containing the drop-seq count matrix (gene x sample).

	Examples
	--------
	>>> tools.load_dropseq_file('example.umi.dge.txt.gz', 'GRCh38')
	"""

	df = pd.read_csv(input_file, header = 0, index_col = 0, sep = '\t', compression = 'gzip')
	channel = {"barcodes" : np.array([x.encode() for x in df.columns]), 
			   "matrix" : df}
	results = {genome : channel}

	return results



def write_10x_h5_file(output_h5, genome2data):
	"""Write count matricies into one 10x v2 formatted hdf5 file output_h5

	Parameters
	----------

	output_h5 : `str`
		The output file name.
	genome2data : `dict of dict`
		This is a dictionary containing genome-channel pairs. Each channel is a dictionary containing one count matrix.

	Returns
	-------

	None

	Examples
	--------
	>>> tools.write_10x_h5_file('example_10x.h5', results_map)
	"""

	with tables.open_file(output_h5, mode="w", title=output_h5, filters=tables.Filters(complevel=1)) as hd5_out:
		for genome, output_data in genome2data.items():
			out_group = hd5_out.create_group("/", genome)
			hd5_out.create_carray(out_group, "barcodes", obj=output_data["barcodes"])
			hd5_out.create_carray(out_group, "data", obj=output_data["matrix"].data)
			hd5_out.create_carray(out_group, "indices", obj=output_data["matrix"].indices)
			hd5_out.create_carray(out_group, "indptr", obj=output_data["matrix"].indptr)
			M, N = output_data["matrix"].shape
			hd5_out.create_carray(out_group, "shape", obj=(N, M))

			for attr, obj in output_data.items():
				if attr not in excluded:
					hd5_out.create_carray(out_group, attr, obj=obj)



class aggr_matrix:
	def __init__(self, channel, nsample, is_dropseq):
		self.out_hd5 = {}
		self.attrs = set()
		for attr, value in channel.items():
			self.attrs.add(attr)
			self.out_hd5[attr] = value if attr in col_attrs else [value]
		self.nsample = nsample
		self.is_dropseq = is_dropseq is not None

	def append(self, channel, nsample):
		attrs_c = set(channel)
		attrs_r = self.attrs & attrs_c # joint	
		for attr in attrs_r:
			if attr not in col_attrs:
				self.out_hd5[attr].append(channel[attr])
		attrs_r = self.attrs - attrs_c # only in existing
		for attr in attrs_r:
			self.out_hd5[attr].append(np.repeat('', nsample))
		attrs_r = attrs_c - self.attrs # only in new
		for attr in attrs_r:
			self.out_hd5[attr] = [np.repeat('', self.nsample), channel[attr]]
		self.nsample += nsample
		self.attrs = self.attrs | attrs_c

	def merge(self):
		if not self.is_dropseq:
			self.out_hd5["matrix"] = vstack(self.out_hd5["matrix"], "csr")
		else:
			df_new = pd.concat(self.out_hd5["matrix"], axis = 1, sort = True)
			df_new.fillna(0, inplace = True)
			self.out_hd5["matrix"] = csr_matrix(df_new.astype(int).T.values)
			self.out_hd5["genes"] = self.out_hd5["gene_names"] = df_new.index.values.astype(str)

		self.out_hd5["barcodes"] = np.concatenate(self.out_hd5["barcodes"])

		for attr in self.attrs:
			if (attr not in col_attrs) and (attr not in excluded):
				self.out_hd5[attr] = np.concatenate(self.out_hd5[attr])



def find_digits(value):
	pos = len(value) - 1
	while pos >= 0 and value[pos].isdigit():
		pos -= 1
	pos += 1
	assert pos < len(value)
	return (value[:pos], int(value[pos:]))


def parse_restriction_string(rstr):
	pos = rstr.index(':')
	name = rstr[: pos]
	isin = True
	if rstr[pos + 1] == '~':
		isin = False
		pos += 1
	content = set()
	for item in rstr[pos + 1:].split(','):
		values = item.split('-')
		if len(values) == 1:
			content.add(values[0])
		else:
			prefix, fr = find_digits(values[0])
			assert values[1].isdigit()
			to = int(values[1]) + 1
			for i in range(fr, to):
				content.add(prefix + str(i))
	return (name, isin, content)

def aggregate_10x_matrices(csv_file, restrictions, attributes, output_file, google_cloud = False, select_singlets = False, ngene = None, is_dropseq = None):
	"""Aggregate channel-specific 10x count matrices into one big count matrix.

	This function takes as input a csv_file, which contains at least 2 columns — Sample, sample name; Location, folder that contains the count matrices (e.g. filtered_gene_bc_matrices_h5.h5). It outputs a 10x-formatted HDF5 file containing on merged matrix per genome.
	
	Parameters
	----------

	csv_file : `str`
		The CSV file containing information about each 10x channel.
	restrictions : `list[str]`
		A list of restrictions used to select channels, each restriction takes the format of name:value,…,value or name:~value,..,value, where ~ refers to not.
	attributes : `list[str]`
		A list of attributes need to be incorporated into the output count matrix.
	output_file : `str`
		The output count matrix file name.
	google_cloud : `bool`, optional (default: False) 
		If the channel-specific count matrices are stored in a google bucket.
	select_singlets : `bool`, optional (default: False)
		If we have demultiplexed data, turning on this option will make scCloud only include barcodes that are predicted as singlets.
	ngene : `int`, optional (default: None)
		The minimum number of expressed genes to keep one barcode. 
	is_dropseq: `str`, optional (default: None)
		If None, assume data are generated from 10 assays; otherwise, assume data are generated using Drop-Seq with reference genome <is_dropseq>.

	Returns
	-------
	
	None

	Examples
	--------
	>>> tools.aggregate_matrix('example.csv', ['Source:pbmc', 'Donor:1'], ['Source', 'Platform', 'Donor'], 'example_10x.h5')
	"""

	df = pd.read_csv(csv_file, header=0, index_col='Sample')
	df['Sample'] = df.index

	# Automatically detect assay type
	input_type = 'dropseq' if df['Location'].iat[0].endswith('dge.txt.gz') else '10x' 

	# Select channels
	rvec = [parse_restriction_string(x) for x in restrictions]
	
	if attributes is None:
		attributes = []
	
	idx = pd.Series([True] * df.shape[0], index=df.index, name='selected')
	for name, isin, content in rvec:
		assert name in df.columns
		if isin:
			idx = idx & df[name].isin(content)
		else:
			idx = idx & (~(df[name].isin(content)))

	df = df.loc[idx]

	if df.shape[0] == 0:
		print("No channels pass the restrictions!")
		return

	# Load channels
	tot = 0
	merged_map = {}
	for sample_name, row in df.iterrows():
		input_hd5_file = row['Location']
		if google_cloud:
			call_args = ['gsutil', '-m', 'cp', input_hd5_file, '{0}_tmp_h5.h5'.format(sample_name)]
			check_call(call_args)
			input_hd5_file = '{0}_tmp_h5.h5'.format(sample_name)
		
		results = load_10x_h5_file(input_hd5_file, ngene, select_singlets) if is_dropseq is None else load_dropseq_file(input_hd5_file, is_dropseq)
		for genome, channel in results.items():
			barcodes = np.array(['{}-{}'.format(sample_name, y[:-2] if y.endswith("-1") else y) for y in (x.decode() for x in channel["barcodes"])])
			nsample = barcodes.size
			channel["barcodes"] = barcodes
			for attr in attributes:
				channel[attr] = np.repeat(row[attr], nsample)
			if genome in merged_map:
				merged_map[genome].append(channel, nsample)
			else:
				merged_map[genome] = aggr_matrix(channel, nsample, is_dropseq)

		tot += 1
		print("Processed {}.".format(input_hd5_file))

	# Delete temporary file
	if google_cloud: 
		check_call(['rm', '-f', input_hd5_file])

	# Merge channels
	results_map = {}
	for genome, aggr in merged_map.items():
		aggr.merge()
		results_map[genome] = aggr.out_hd5

	write_10x_h5_file(output_file, results_map)
	print("Generated {file} from {tot} files.".format(file=output_file, tot=tot))
