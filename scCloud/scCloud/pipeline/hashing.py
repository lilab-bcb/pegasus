import os
import pandas as pd
from .. import tools, demuxEM

def run_hashing_pipeline(input_adt_file, input_rna_file, output_name, **kwargs):
	# determine genome if necessary
	if kwargs['genome'] is None:
		kwargs['genome'] = tools.determine_genome(input_rna_file)

	# load input data
	adt = tools.read_input(input_adt_file)
	data = tools.read_input(input_rna_file, genome = kwargs['genome'], ngene = kwargs['min_num_genes'])

	# run demuxEM
	demuxEM.demultiplex(data, adt, unknown = kwargs['unknown_rate'], alpha_value = kwargs['alpha_value'], n_threads = kwargs['n_jobs'])

	# annotate raw matrix with demuxEM results
	raw_data = tools.load_10x_h5_file(input_rna_file, kwargs['genome'], threshold = None)

	barcodes = pd.Index(raw_data['barcodes'])
	idx = barcodes.isin(data.obs_names)
	selected = barcodes[idx]

	demux_type = np.empty(barcodes.size, dtype = 'object')
	demux_type[:] = ''
	demux_type[idx] = data.obs.loc[selected, 'demux_type']
	raw_data['demux_type'] = demux_type.astype('str')

	assignment = np.empty(barcodes.size, dtype = 'object')
	assignment[:] = ''
	assignment[idx] = data.obs.loc[selected, 'assignment']
	raw_data['assignment'] = assignment.astype('str')

	# output results
	tools.write_10x_h5_file(output_name + '_demux_10x.h5', raw_data, kwargs['genome'], ['demux_type', 'assignment'])
	adata.write(output_name + "_demux.h5ad")
	print("Results are written.")
