import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scCloud.tools import load_10x_h5_file, write_10x_h5_file



def load_antibody_csv(input_csv, antibody_control_csv):
	barcodes = []
	antibody_names = []
	antibody_to_pos = {}
	stacks = []

	with open(input_csv) as fin:
		barcodes = np.array([(x + '-1').encode() for x in next(fin).strip().split(',')[1:]])
		for i, line in enumerate(fin):
			fields = line.strip().split(',')
			antibody_names.append(fields[0])
			antibody_to_pos[fields[0]] = i
			stacks.append([int(x) for x in fields[1:]])

	adt_matrix = np.stack(stacks, axis = 1).astype(float)
	antibody_names = np.array([x.encode() for x in antibody_names])

	if antibody_control_csv is None:
		adt_matrix = np.log(adt_matrix + 1.0)
	else:
		series = pd.read_csv(antibody_control_csv, header = 0, index_col = 0, squeeze = True)
		idx = np.zeros(len(antibody_names), dtype = bool)
		for antibody, control in series.iteritems():
			pos_a = antibody_to_pos[antibody]
			pos_c = antibody_to_pos[control]
			idx[pos_a] = True
			# convert to log expression
			adt_matrix[:, pos_a] = np.maximum(np.log(adt_matrix[:, pos_a] + 1.0) - np.log(adt_matrix[:, pos_c] + 1.0), 0.0)
		adt_matrix = adt_matrix[:, idx]
		antibody_names = antibody_names[idx]

	channel = {"barcodes" : barcodes, "antibody_names" : antibody_names, "matrix" : csr_matrix(adt_matrix)}
	return channel

def merge_rna_and_adt_data(input_raw_h5, input_csv, antibody_control_csv, output_10x_h5):
	results = load_10x_h5_file(input_raw_h5)
	print("Loaded the RNA matrix.")
	assert len(results) == 1
	genome = next(iter(results))
	results['CITE_Seq_' + genome] = load_antibody_csv(input_csv, antibody_control_csv)
	print("Loaded the ADT matrix.")
	write_10x_h5_file(output_10x_h5, results)
	print("Merged output is written.")

def capping(adt_matrix, percentile):
	for i in range(adt_matrix.shape[1]):
		cap = np.percentile(adt_matrix[:, i], percentile)
		adt_matrix[adt_matrix[:, i] > cap, i] = cap
