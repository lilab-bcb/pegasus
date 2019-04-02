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

from . import calculate_affinity_matrix, calculate_nearest_neighbors, select_cells, construct_graph



def calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state, init = 'random', n_iter = 1000, n_iter_early_exag = 250):
	tsne = TSNE(n_jobs = n_jobs, n_components = n_components, perplexity = perplexity, early_exaggeration = early_exaggeration, learning_rate = learning_rate, random_state = random_state, verbose = 1, init = init, n_iter = n_iter, n_iter_early_exag = n_iter_early_exag)
	X_tsne = tsne.fit_transform(X)
	print("Final error = {}".format(tsne.kl_divergence_))
	return X_tsne

def calc_fitsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state):
	# FItSNE will change X content
	return FItSNE(X.astype('float64'), nthreads = n_jobs, no_dims = n_components, perplexity = perplexity, early_exag_coeff = early_exaggeration, learning_rate = learning_rate, rand_seed = random_state)

def calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state):
	umap = UMAP(n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, spread = spread, random_state = random_state)
	return umap.fit_transform(X)

def calc_force_directed_layout(W, file_name, n_jobs, target_change_per_node, target_steps, is3d, memory, random_state):
	input_graph_file = '{file_name}.net'.format(file_name = file_name)
	G = construct_graph(W, adjust_weights = True)
	G.write(input_graph_file)
	print(input_graph_file + ' is written.')			

	output_coord_file = '{file_name}.coords.txt'.format(file_name = file_name)

	classpath = pkg_resources.resource_filename('scCloud', 'ext/forceatlas2.jar') + ':' + pkg_resources.resource_filename('scCloud', 'ext/gephi-toolkit-0.9.2-all.jar')
	command = ['java', '-Djava.awt.headless=true', '-Xmx{memory}g'.format(memory = memory), '-cp', classpath, \
				'kco.forceatlas2.Main', '--input', input_graph_file, '--output', file_name + '.coords', \
				'--nthreads', str(n_jobs), '--seed', str(random_state), '--targetChangePerNode', str(target_change_per_node), '--targetSteps', str(target_steps)]
	if not is3d:
		command.append('--2d')
	check_call(command)
	print("Force-directed layout is generated.")

	fle_coords = pd.read_csv(output_coord_file, header = 0, index_col = 0, sep = '\t').values
	check_call(['rm', '-f', input_graph_file])
	check_call(['rm', '-f', output_coord_file])

	return fle_coords



def run_tsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 100, out_basis = 'tsne'):
	start = time.time()
	X_tsne = calc_tsne(data.obsm[rep_key], n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	data.obsm['X_' + out_basis] = X_tsne
	end = time.time()
	print("tSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 100, out_basis = 'fitsne'):
	start = time.time()
	data.obsm['X_' + out_basis] = calc_fitsne(data.obsm[rep_key], n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state) 
	end = time.time()
	print("FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 100, out_basis = 'umap'):
	start = time.time()
	data.obsm['X_' + out_basis] = calc_umap(data.obsm[rep_key], n_components, n_neighbors, min_dist, spread, random_state)
	end = time.time()
	print("UMAP is calculated. Time spent = {:.2f}s.".format(end - start))

def run_force_directed_layout(data, file_name, n_jobs, K = 50, target_change_per_node = 2.0, target_steps = 10000, is3d = False, memory = 8, random_state = 100, out_basis = 'fle'):
	start = time.time()
	K = min(K - 1, data.uns['diffmap_knn_indices'].shape[1]) # K - 1: exclude self
	W = calculate_affinity_matrix(data.uns['diffmap_knn_indices'][:, 0:K], data.uns['diffmap_knn_distances'][:, 0:K])
	data.obsm['X_' + out_basis] = calc_force_directed_layout(W, file_name, n_jobs, target_change_per_node, target_steps, is3d, memory, random_state);
	end = time.time()
	print("Force-directed layout is calculated. Time spent = {:.2f}s.".format(end - start))



def run_net_tsne(data, rep_key, selected, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 100, net_alpha = 0.1, polish_learning_rate = 1e5, polish_n_iter = 150, out_basis = 'net_tsne'):
	start = time.time()
	X = data.obsm[rep_key][selected,:].astype('float64')
	X_tsne = calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = net_alpha, random_state = random_state)
	regressor.fit(X, X_tsne)
	X_tsne_pred = regressor.predict(data.obsm[rep_key][~selected,:].astype('float64'))
	Y_init = np.zeros((data.shape[0], 2), dtype = np.float64)
	Y_init[selected,:] = X_tsne
	Y_init[~selected,:] = X_tsne_pred
	data.obsm['X_' + out_basis] = calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, polish_learning_rate, random_state, init = Y_init, n_iter = polish_n_iter, n_iter_early_exag = 0)
	end = time.time()
	print("Net tSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 100, knn_indices = 'diffmap_knn_indices', first_K = 5):
	start = time.time()
	selected = select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	X = data.obsm[rep_key][selected,:].astype('float64')
	X_fitsne = calc_fitsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01, random_state = random_state)
	regressor.fit(X, X_fitsne)
	X_fitsne_pred = regressor.predict(data.obsm[rep_key][~selected,:].astype('float64'))
	data.obsm['X_net_fitsne'] = np.zeros((data.shape[0], 2), dtype = np.float64)
	data.obsm['X_net_fitsne'][selected,:] = X_fitsne
	data.obsm['X_net_fitsne'][~selected,:] = X_fitsne_pred
	end = time.time()
	print("Net FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 100, knn_indices = 'diffmap_knn_indices', first_K = 5):
	start = time.time()
	selected = select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	X = data.obsm[rep_key][selected,:].astype('float64')
	X_umap = calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state)
	regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01, random_state = random_state)
	regressor.fit(X, X_umap)
	X_umap_pred = regressor.predict(data.obsm[rep_key][~selected,:].astype('float64'))
	data.obsm['X_net_umap'] = np.zeros((data.shape[0], 2), dtype = np.float64)
	data.obsm['X_net_umap'][selected,:] = X_umap
	data.obsm['X_net_umap'][~selected,:] = X_umap_pred
	end = time.time()
	print("Net UMAP is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_fle(data, file_name, n_jobs, K = 50, layout = 'fa', n_steps = 10000, memory = 20, random_state = 100, knn_indices = 'diffmap_knn_indices', first_K = 5):
	# start = time.time()
	# selected = select_cells(data.uns[knn_indices], first_K, random_state = random_state)
	# X = data.obsm['X_diffmap'][selected,:]
	# indices, distances, knn_index = calculate_nearest_neighbors(X, n_jobs, K = K, random_state = random_state, full_speed = False)
	# W = calculate_affinity_matrix(indices, distances)
	# X_fle = calc_force_directed_layout(W, file_name, n_jobs, layout, n_steps, memory)
	# regressor = MLPRegressor(hidden_layer_sizes = (100, 70, 50, 25), activation = 'relu', solver = 'sgd', learning_rate = 'adaptive', alpha = 0.01, random_state = random_state)
	# regressor.fit(X, X_fle)
	# X_fle_pred = regressor.predict(data.obsm['X_diffmap'][~selected,:])
	# data.obsm['X_net_fle'] = np.zeros((data.shape[0], 2))
	# data.obsm['X_net_fle'][selected,:] = X_fle
	# data.obsm['X_net_fle'][~selected,:] = X_fle_pred
	# end = time.time()
	print("Net FLE is calculated. Time spent = {:.2f}s.".format(end - start))


