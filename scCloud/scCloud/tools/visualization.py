import time
import numpy as np
import pandas as pd
from subprocess import check_call

import os
import pkg_resources

from MulticoreTSNE import MulticoreTSNE as TSNE
from umap import UMAP
from fitsne import FItSNE

def run_tsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, out_basis = 'tsne'):
	start = time.time()
	tsne = TSNE(n_jobs = n_jobs, n_components = n_components, perplexity = perplexity, early_exaggeration = early_exaggeration, learning_rate = learning_rate, random_state = random_state)
	X = data.obsm[rep_key].astype('float64')
	X_tsne = tsne.fit_transform(X)
	data.obsm['X_' + out_basis] = X_tsne
	end = time.time()
	print("TSNE is done. Time spent = {:.2f}s.".format(end - start))

def run_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, random_state = 0):
	start = time.time()
	X = data.obsm[rep_key].astype('float64').copy(order = 'C')
	X_fitsne = FItSNE(X, nthreads = n_jobs, no_dims = n_components, perplexity = perplexity, rand_seed = (random_state if random_state is not None else -1), early_exag_coeff = early_exaggeration)
	data.obsm['X_fitsne'] = X_fitsne
	end = time.time()
	print("FItSNE is done. Time spent = {:.2f}s.".format(end - start))

def run_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 0):
	start = time.time()
	umap = UMAP(n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, spread = spread, random_state = random_state)
	X = data.obsm[rep_key].astype('float64')
	X_umap = umap.fit_transform(X)
	data.obsm['X_umap'] = X_umap	
	end = time.time()
	print("UMAP is done. Time spent = {:.2f}s.".format(end - start))

def run_force_directed_layout(data, file_name, n_jobs, affinity = 'W_diffmap', K = 50, layout = 'fa', n_steps = 10000, memory = 20):
	start = time.time()
	W = data.uns[affinity]

	input_graph_file = '{file_name}.net'.format(file_name = file_name)
	output_coord_file = '{file_name}.coords.txt'.format(file_name = file_name)

	with open(input_graph_file, 'w') as writer:	
		n_obs = W.shape[0]
		writer.write("*Vertices {n_obs}\n".format(n_obs = n_obs))
		for i in range(n_obs):
			writer.write("{node} \"{node}\"\n".format(node = i + 1))
		writer.write("*Edges\n")
		rows, cols = W.nonzero()
		for i, j in zip(rows, cols):
			if i < j:
				writer.write("{u} {v} {w:.6g}\n".format(u = i + 1, v = j + 1, w = W[i, j]))
	print(input_graph_file + ' is written.')			

	classpath = os.path.dirname(pkg_resources.resource_filename('scCloud', 'ext/GraphLayout.class')) + ':' + \
				pkg_resources.resource_filename('scCloud', 'ext/gephi-toolkit-0.9.2-all.jar')
	check_call(['java', '-Djava.awt.headless=true', '-Xmx{memory}g'.format(memory = memory), '-cp', classpath, \
				'GraphLayout', input_graph_file, output_coord_file, layout, str(n_steps), str(n_jobs)])
	print("Force-directed layout is generated.")

	data.obsm['X_fle'] = pd.read_table(output_coord_file, header = 0, index_col = 0).values
	print("X_fle is written.")

	check_call(['rm', '-f', input_graph_file])
	check_call(['rm', '-f', output_coord_file])

	end = time.time()
	print("Force-directed layout is done. Time spent = {:.2f}s.".format(end - start))



