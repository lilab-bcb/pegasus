import time
import numpy as np
import pandas as pd
from subprocess import check_call

import os
import pkg_resources

from MulticoreTSNE import MulticoreTSNE as TSNE
from umap import UMAP
from fitsne import FItSNE
from sklearn.neural_network import MLPRegressor

from . import calculate_affinity_matrix, calculate_nearest_neighbors, select_cells


def calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state):
	tsne = TSNE(n_jobs = n_jobs, n_components = n_components, perplexity = perplexity, early_exaggeration = early_exaggeration, learning_rate = learning_rate, random_state = random_state)
	return tsne.fit_transform(X)

def calc_fitsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state):
	return FItSNE(X, nthreads = n_jobs, no_dims = n_components, perplexity = perplexity, early_exag_coeff = early_exaggeration, learning_rate = learning_rate, rand_seed = (random_state if random_state is not None else -1))

def calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state):
	umap = UMAP(n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, spread = spread, random_state = random_state)
	return umap.fit_transform(X)

def calc_force_directed_layout(W, file_name, n_jobs, n_steps, memory):
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

	fle_coords = pd.read_table(output_coord_file, header = 0, index_col = 0).values
	check_call(['rm', '-f', input_graph_file])
	check_call(['rm', '-f', output_coord_file])

	return fle_coords



def run_tsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0):
	start = time.time()
	X = data.obsm[rep_key].astype('float64')
	X_tsne = calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	data.obsm['X_tsne'] = X_tsne
	end = time.time()
	print("tSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0):
	start = time.time()
	X = data.obsm[rep_key].astype('float64').copy(order = 'C')
	data.obsm['X_fitsne'] = calc_fitsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	end = time.time()
	print("FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 0):
	start = time.time()
	X = data.obsm[rep_key].astype('float64')
	data.obsm['X_umap'] = calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state)
	end = time.time()
	print("UMAP is calculated. Time spent = {:.2f}s.".format(end - start))

def run_force_directed_layout(data, file_name, n_jobs, K = 50, n_steps = 10000, memory = 20):
	start = time.time()
	K = min(K - 1, data.uns['diffmap_knn_indices'].shape[1]) # K - 1: exclude self
	W = calculate_affinity_matrix(data.uns['diffmap_knn_indices'][:, 0:K], data.uns['diffmap_knn_distances'][:, 0:K])
	data.obsm['X_fle'] = calc_force_directed_layout(W, file_name, n_jobs, n_steps, memory)
	end = time.time()
	print("Force-directed layout is calculated. Time spent = {:.2f}s.".format(end - start))



def run_net_tsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, knn_indices = 'diffmap_knn_indices', first_K = 5):
	start = time.time()
	selected = scCloud.tools.select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	X = data.obsm[rep_key][selected,:].astype('float64')
	X_tsne = calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01)
	regressor.fit(X, X_tsne)
	X_tsne_pred = regressor.predict(data.obsm[rep_key][~selected,:].astype('float64'))
	data.obsm['X_net_tsne'] = np.zeros((data.shape[0], 2), dtype = np.float64)
	data.obsm['X_net_tsne'][selected,:] = X_tsne
	data.obsm['X_net_tsne'][~selected,:] = X_tsne_pred
	end = time.time()
	print("Net tSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, knn_indices = 'diffmap_knn_indices', first_K = 5):
	start = time.time()
	selected = scCloud.tools.select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	X = data.obsm[rep_key][selected,:].astype('float64').copy(order = 'C')
	X_fitsne = calc_fitsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01)
	regressor.fit(X, X_fitsne)
	X_fitsne_pred = regressor.predict(data.obsm[rep_key][~selected,:].astype('float64'))
	data.obsm['X_net_fitsne'] = np.zeros((data.shape[0], 2), dtype = np.float64)
	data.obsm['X_net_fitsne'][selected,:] = X_fitsne
	data.obsm['X_net_fitsne'][~selected,:] = X_fitsne_pred
	end = time.time()
	print("Net FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 0, knn_indices = 'diffmap_knn_indices', first_K = 5):
	start = time.time()
	selected = scCloud.tools.select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	X = data.obsm[rep_key][selected,:].astype('float64')
	X_umap = calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01)
	regressor.fit(X, X_umap)
	X_umap_pred = regressor.predict(data.obsm[rep_key][~selected,:].astype('float64'))
	data.obsm['X_net_umap'] = np.zeros((data.shape[0], 2), dtype = np.float64)
	data.obsm['X_net_umap'][selected,:] = X_umap
	data.obsm['X_net_umap'][~selected,:] = X_umap_pred
	end = time.time()
	print("Net UMAP is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_fle(data, file_name, n_jobs, K = 50, n_steps = 10000, memory = 20, random_state = 0, knn_indices = 'diffmap_knn_indices', first_K = 5):
	start = time.time()
	selected = scCloud.tools.select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	X = data.obsm['X_diffmap'][selected,:]
	indices, distances, knn_index = calculate_nearest_neighbors(X, n_jobs, K = K, random_state = random_state, full_speed = False)
	W = calculate_affinity_matrix(indices, distances)
	X_fle = calc_force_directed_layout(W, file_name, n_jobs, n_steps, memory)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01)
	regressor.fit(X, X_fle)
	X_fle_pred = regressor.predict(data.obsm['X_diffmap'][~selected,:])
	data.obsm['X_net_fle'] = np.zeros((data.shape[0], 2))
	data.obsm['X_net_fle'][selected,:] = X_fle
	data.obsm['X_net_fle'][~selected,:] = X_fle_pred
	end = time.time()
	print("Net FLE is calculated. Time spent = {:.2f}s.".format(end - start))
