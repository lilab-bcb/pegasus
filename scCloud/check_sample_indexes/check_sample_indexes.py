#!/usr/bin/env python

from sys import exit

import json
import pkg_resources

def hamming(index1, index2):
	return sum([index1[i] != index2[i] for i in range(min(len(index1), len(index2)))])

def validate_indices(index_arr, n_mis = 1):
	min_hamming = n_mis * 2 + 1

	s = len(index_arr)
	pass_validation = True
	for i in range(s - 1):
		for j in range(i + 1, s):
			hd = hamming(index_arr[i], index_arr[j])
			if hd < min_hamming:
				print("Warning: index {0} and {1} have a hamming distance of {2} < {3}!".format(index_arr[i], index_arr[j], hd, min_hamming))
				pass_validation = False

	return pass_validation

def check_index_set(index_set, index_arr, n_mis = 1):
	min_hamming = n_mis * 2 + 1

	for sequence in index_set[1]:
		for idx in index_arr:
			hd = hamming(sequence, idx)
			if hd < min_hamming:
				return False

	return True

def load_index_file(index_file):
	# Load index file
	index_arr = []
	with open(index_file) as fin:
		for line in fin:
			index_arr.append(line.strip())
	return index_arr

def load_chromium_indexes():
	# Load chromium index sets
	sample_indexes_file = pkg_resources.resource_filename('scCloud.check_sample_indexes', 'chromium-dna-sample-indexes-plate.json')
	with open(sample_indexes_file) as fin:
		return json.load(fin)

def run_check_sample_indexes(index_file, n_mis = 1, n_report = 9999):
	index_arr = load_index_file(index_file)

	if not validate_indices(index_arr, n_mis):
		print("Error: Index collision detected in {0}!".format(index_file))
		exit(-1)

	indexes = load_chromium_indexes()

	n_success = 0
	for index_set in indexes:
		if check_index_set(index_set, index_arr, n_mis):
			print(index_set[0])
			n_success += 1
			if n_success >= n_report:
				break
