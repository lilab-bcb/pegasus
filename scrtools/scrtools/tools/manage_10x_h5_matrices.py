#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from scipy.sparse import csc_matrix, hstack
import tables
import copy
from subprocess import check_call


basic_attrs = set(["barcodes", "gene_names", "genes", "matrix"])

def load_10x_h5_file(input_h5, genome):
	inpmat = {}
	with tables.open_file(input_h5) as h5_in:
		for node in h5_in.walk_nodes("/" + genome, "Array"):
			inpmat[node.name] = node.read()

		inpmat["matrix"] = csc_matrix((inpmat["data"], inpmat["indices"], inpmat["indptr"]), shape=inpmat["shape"])

		inpmat.pop("data")
		inpmat.pop("indices")
		inpmat.pop("indptr")
		inpmat.pop("shape")

	attributes = [key for key in inpmat.keys() if key not in basic_attrs]

	return inpmat, attributes


def write_10x_h5_file(output_h5, output_data, genome, attributes):
	with tables.open_file(output_h5, mode="w", title=output_h5, filters=tables.Filters(complevel=1)) as hd5_out:
		out_group = hd5_out.create_group("/", genome)
		hd5_out.create_carray(out_group, "barcodes", obj=output_data["barcodes"])
		hd5_out.create_carray(out_group, "gene_names", obj=output_data["gene_names"])
		hd5_out.create_carray(out_group, "genes", obj=output_data["genes"])
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


def aggregate_10x_matrices(csv_file, genome, restrictions, attributes, output_file, google_cloud=False):
	df = pd.read_csv(csv_file, header=0, index_col='Sample')
	restrictions.append("Reference:{}".format(genome))
	rvec = [parse_restriction_string(x) for x in restrictions]
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
			call_args = ['gsutil', '-m', 'cp', input_hd5_file,
						 '{0}_filtered_gene_bc_matrices_h5.h5'.format(sample_name)]
			check_call(call_args)
			input_hd5_file = '{0}_filtered_gene_bc_matrices_h5.h5'.format(sample_name)
		
		channel, tmp_attrs = load_10x_h5_file(input_hd5_file, genome)
		channel["barcodes"] = np.array([(sample_name + '-' + x.decode()).encode() for x in channel["barcodes"]])

		if out_hd5 is None:
			out_hd5 = copy.copy(channel)
			out_hd5["matrix"] = []
			out_hd5["barcodes"] = []
			for attr in attributes:
				out_hd5[attr] = []

		out_hd5["matrix"].append(channel["matrix"])
		out_hd5["barcodes"].append(channel["barcodes"])
		for attr in attributes:
			if attr in channel:
				out_hd5[attr].append(channel[attr])
			else:
				out_hd5[attr].append(np.repeat(row[attr], channel["matrix"].shape[1]))

		print("Processed {}.".format(input_hd5_file))
		tot += 1

	out_hd5["matrix"] = hstack(out_hd5["matrix"], "csc")
	out_hd5["barcodes"] = np.concatenate(out_hd5["barcodes"])
	for attr in attributes:
		out_hd5[attr] = np.concatenate(out_hd5[attr])

	write_10x_h5_file(output_file, out_hd5, genome, attributes)
	print("Generated {file} from {tot} files.".format(file=output_file, tot=tot))
