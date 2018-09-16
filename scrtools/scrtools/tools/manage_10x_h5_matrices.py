#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from scipy.sparse import csc_matrix, hstack
import tables
import copy
from subprocess import check_call



def load_10x_h5_file(input_h5, genome, threshold = 30000, ngene = 100):
	"""Load 10x format matrix from h5 file
	
	Parameters
	----------

	input_h5 : `str`
		The matrix in h5 format.
	genome : `str`
		Genome string.
	threshold : `int`, optional (default: 30000)
		If matrix contain more than threshold barcodes, filter barcodes with low number of genes.
	ngene : `int`, optional (default: 100)
		Minimum number of genes to keep a barcode if # of barcodes > threshold.
	"""	
	inpmat = {}
	with tables.open_file(input_h5) as h5_in:
		for node in h5_in.walk_nodes("/" + genome, "Array"):
			inpmat[node.name] = node.read()

	mat = csc_matrix((inpmat["data"], inpmat["indices"], inpmat["indptr"]), shape=inpmat["shape"])
	if mat.shape[1] > threshold:
		selected = mat.getnnz(axis = 0) >= ngene
		mat = mat[:, selected]
		inpmat["barcodes"] = inpmat["barcodes"][selected]

	inpmat["matrix"] = mat

	inpmat.pop("data")
	inpmat.pop("indices")
	inpmat.pop("indptr")
	inpmat.pop("shape")

	return inpmat


def load_antibody_csv(input_csv):
	barcodes = []
	antibody_names = []
	stacks = []
	with open(input_csv) as fin:
		barcodes = next(fin).strip().split(',')[1:]
		for line in fin:
			fields = line.strip().split(',')
			antibody_names.append(fields[0])
			stacks.append([int(x) for x in fields[1:]])

	inpmat = {}
	inpmat["barcodes"] = [x.encode() for x in barcodes]
	inpmat["antibody_names"] = [x.encode() for x in antibody_names]
	inpmat["matrix"] = csc_matrix(np.stack(stacks, axis = 0))

	return inpmat


def write_10x_h5_file(output_h5, output_data, genome, attributes):
	with tables.open_file(output_h5, mode="w", title=output_h5, filters=tables.Filters(complevel=1)) as hd5_out:
		out_group = hd5_out.create_group("/", genome)
		hd5_out.create_carray(out_group, "barcodes", obj=output_data["barcodes"])
		hd5_out.create_carray(out_group, "data", obj=output_data["matrix"].data)
		hd5_out.create_carray(out_group, "indices", obj=output_data["matrix"].indices)
		hd5_out.create_carray(out_group, "indptr", obj=output_data["matrix"].indptr)
		hd5_out.create_carray(out_group, "shape", obj=output_data["matrix"].shape)

		for attr in attributes:
			hd5_out.create_carray(out_group, attr, obj=output_data[attr])


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


def aggregate_10x_matrices(csv_file, genome, restrictions, attributes, output_name, google_cloud=False, input_type='gene'):
	"""Aggregate channel-specific 10x count matrices into one big count matrix.

	This function takes as input a csv_file, which contains at least 3 columns — Sample, sample name; Location, folder that contains the count matrices (e.g. filtered_gene_bc_matrices_h5.h5); Reference, genome reference used for 10x cellranger. It outputs a 10x-formatted HDF5 file for the big count matrix.
	
	Parameters
	----------

	csv_file : `str`
		The CSV file containing information about each 10x channel.
	genome : `str`
		The genome each sample comes from.
	restrictions : `list[str]`
		A list of restrictions used to select channels, each restriction takes the format of name:value,…,value or name:~value,..,value, where ~ refers to not.
	attributes : `list[str]`
		A list of attributes need to be incorporated into the output count matrix.
	output_name : `str`
		The output count matrix file name prefix. If input_type == 'gene', output_name_10x.h5 will be generated. If input_type == 'ADT', output_name.h5at will be generated.
	google_cloud : `bool`, optional (default: `False`) 
		If the channel-specific count matrices are stored in a google bucket.
	input_type : `str`, optional (default: `gene`)
		Input type, 'gene' refers to 10x h5 format; 'dropseq' refers to drop-seq format; 'ADT' refers to CITE-Seq csv.

	Returns
	-------
	
	None

	Examples
	--------
	>>> tools.aggregate_matrix('example.csv', 'GRCh38', ['Source:pbmc', 'Donor:1'], ['Source', 'Platform', 'Donor'], 'example_10x.h5')
	"""

	df = pd.read_csv(csv_file, header=0, index_col='Sample')
	df['Sample'] = df.index
	rvec = [parse_restriction_string(x) for x in restrictions]
	rvec.append(('Reference', True, {genome}))
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

	tot = 0
	out_hd5 = None
	for sample_name, row in df.iterrows():
		input_hd5_file = row['Location']
		if google_cloud:
			call_args = ['gsutil', '-m', 'cp', input_hd5_file, '{0}_tmp_h5.h5'.format(sample_name)]
			check_call(call_args)
			input_hd5_file = '{0}_tmp_h5.h5'.format(sample_name)
		
		if input_type == 'dropseq':
			channel = pd.read_table(input_hd5_file, header = 0, index_col = 0, compression = 'gzip')
			channel.columns = [sample_name + '-' + x for x in channel.columns]

			if out_hd5 is None:
				out_hd5 = {"matrix" : [], "barcodes" : None, "genes" : None, "gene_names" : None}
				for attr in attributes:
					out_hd5[attr] = []

			out_hd5["matrix"].append(channel)
			for attr in attributes:
				out_hd5[attr].append(np.repeat(row[attr], channel.shape[1]))

		else:
			if input_type == 'gene':
				channel = load_10x_h5_file(input_hd5_file, genome)
			elif input_type == 'ADT':
				channel = load_antibody_csv(input_hd5_file)
			else:
				assert(False)

			channel["barcodes"] = np.array([(sample_name + '-' + x.decode()) for x in channel["barcodes"]])

			if out_hd5 is None:
				out_hd5 = copy.copy(channel)
				out_hd5["matrix"] = []
				out_hd5["barcodes"] = []
				for attr in attributes:
					out_hd5[attr] = []

			out_hd5["matrix"].append(channel["matrix"])
			out_hd5["barcodes"].append(channel["barcodes"])
			for attr in attributes:
				out_hd5[attr].append(np.repeat(row[attr], channel["matrix"].shape[1]))

		print("Processed {}.".format(input_hd5_file))
		tot += 1

	# delete temporary file
	if google_cloud: 
		check_call(['rm', '-f', input_hd5_file])

	output_file = output_name
	
	if input_type == 'dropseq':
		df_new = pd.concat(out_hd5["matrix"], axis = 1, sort = True)
		df_new.fillna(0, inplace = True)
		out_hd5["matrix"] = csc_matrix(df_new.astype(int).values)
		out_hd5["barcodes"] = df_new.columns.values.astype(str)
		out_hd5["genes"] = out_hd5["gene_names"] = df_new.index.values.astype(str)
		for attr in attributes:
			out_hd5[attr] = np.concatenate(out_hd5[attr])
		output_file += '_10x.h5'
		attributes.extend(['genes', 'gene_names'])
	else:		
		out_hd5["matrix"] = hstack(out_hd5["matrix"], "csc")
		out_hd5["barcodes"] = np.concatenate(out_hd5["barcodes"])
		for attr in attributes:
			out_hd5[attr] = np.concatenate(out_hd5[attr])
		
		if input_type == 'gene':
			output_file += '_10x.h5'
			attributes.extend(['genes', 'gene_names'])
		elif input_type == 'ADT':
			output_file += '.h5at'
			attributes.append('antibody_names')

	write_10x_h5_file(output_file, out_hd5, genome, attributes)
	print("Generated {file} from {tot} files.".format(file=output_file, tot=tot))
