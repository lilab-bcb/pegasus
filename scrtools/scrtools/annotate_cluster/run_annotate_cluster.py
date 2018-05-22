import anndata
import pkg_resources

from . import annotate_cluster

def run_annotate_cluster(input_file, output_name, threshold, ignoreNA = False, labels = 'louvain_labels'):
	json_file = pkg_resources.resource_filename('scrtools.annotate_cluster', 'cell_type_markers.json')
	
	adata = anndata.read_h5ad(input_file)
	with open(output_name + "_fisher.anno.txt", "w") as fout:
		annotate_cluster.annotate_clusters(adata, 'fisher', json_file, threshold, fout, ignoreNA, labels)
	with open(output_name + "_t.anno.txt", "w") as fout:
		annotate_cluster.annotate_clusters(adata, 't', json_file, threshold, fout, ignoreNA, labels)
