import time
import pkg_resources

from scrtools.tools import read_input
from . import annotate_cluster

def run_annotate_cluster(input_file, output_file, threshold, ignoreNA = False, json_file = "human"):
	start = time.time()
	if json_file == "human_immune":
		json_file = pkg_resources.resource_filename('scrtools.annotate_cluster', 'human_immune_cell_markers.json')
	elif json_file == "mouse_immune":
		json_file = pkg_resources.resource_filename('scrtools.annotate_cluster', 'mouse_immune_cell_markers.json')
	elif json_file == "mouse_brain":
		json_file = pkg_resources.resource_filename('scrtools.annotate_cluster', 'mouse_brain_cell_markers.json')
	adata = read_input(input_file, mode = 'r')
	with open(output_file, 'w') as fout:
		annotate_cluster.annotate_clusters(adata, json_file, threshold, fout, ignoreNA)
	end = time.time()
	print("Time spent for annotating clusters is {:.2f}s.".format(end - start))
