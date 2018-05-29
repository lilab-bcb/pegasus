import time
import anndata
import pkg_resources

from . import annotate_cluster

def run_annotate_cluster(input_file, output_file, threshold, ignoreNA = False):
	start = time.time()
	json_file = pkg_resources.resource_filename('scrtools.annotate_cluster', 'cell_type_markers.json')
	adata = anndata.read_h5ad(input_file, backed = 'r')
	with open(output_file, 'w') as fout:
		annotate_cluster.annotate_clusters(adata, json_file, threshold, fout, ignoreNA)
	end = time.time()
	print("Time spent for annotating clusters is {:.2f}s.".format(end - start))
