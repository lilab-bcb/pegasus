import time
import pkg_resources

from scCloud.tools import read_input
from . import annotate_cluster

def run_annotate_cluster(input_file, output_file, threshold = 0.5, ignoreNA = False, json_file = "human"):
	start = time.time()
	if json_file == "human_immune":
		json_file = pkg_resources.resource_filename('scCloud.annotate_cluster', 'human_immune_cell_markers.json')
	elif json_file == "mouse_immune":
		json_file = pkg_resources.resource_filename('scCloud.annotate_cluster', 'mouse_immune_cell_markers.json')
	elif json_file == "mouse_brain":
		json_file = pkg_resources.resource_filename('scCloud.annotate_cluster', 'mouse_brain_cell_markers.json')
	elif json_file == "human_brain":
		json_file = pkg_resources.resource_filename('scCloud.annotate_cluster', 'human_brain_cell_markers.json')
	data = read_input(input_file, mode = 'r')
	with open(output_file, 'w') as fout:
		annotate_cluster.annotate_clusters(data, json_file, threshold, fout, ignoreNA)
	data.file.close()
	end = time.time()
	print("Time spent for annotating clusters is {:.2f}s.".format(end - start))

def annotate_anndata_object(input_file, annotation):
	data = read_input(input_file, mode = 'r+')
	anno_attr, anno_str = annotation.split(':')
	anno = anno_str.split(';')
	data.obs[anno_attr] = [anno[int(x) - 1] for x in data.obs[data.uns['de_labels']]]
	data.write(input_file)
