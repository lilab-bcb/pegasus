import time
import numpy as np
import pandas as pd
from subprocess import check_call

import os
import pkg_resources

from MulticoreTSNE import MulticoreTSNE as TSNE
from umap import UMAP

from . import calculate_affinity_matrix, calculate_nearest_neighbors, select_cells, construct_graph, net_train_and_predict



def calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state, init = 'random', n_iter = 1000, n_iter_early_exag = 250):
	tsne = TSNE(n_jobs = n_jobs, n_components = n_components, perplexity = perplexity, early_exaggeration = early_exaggeration, learning_rate = learning_rate, \
		random_state = random_state, verbose = 1, init = init, n_iter = n_iter, n_iter_early_exag = n_iter_early_exag)
	X_tsne = tsne.fit_transform(X)
	print("Final error = {}".format(tsne.kl_divergence_))
	return X_tsne


def calc_fitsne(X, nthreads, no_dims, perplexity, early_exag_coeff, learning_rate, rand_seed, initialization = None, max_iter = 1000, stop_early_exag_iter = 250, mom_switch_iter = 250):
	# FItSNE will change X content

	# Check if fftw3 is installed.
	import ctypes.util
	fftw3_loc = ctypes.util.find_library('fftw3')
	if fftw3_loc is None:
		raise Exception("Please install 'fftw3' first to use the FIt-SNE feature!")

	from fitsne import FItSNE

	return FItSNE(X.astype('float64'), nthreads = nthreads, no_dims = no_dims, perplexity = perplexity, early_exag_coeff = early_exag_coeff, learning_rate = learning_rate, \
		rand_seed = rand_seed, initialization = initialization, max_iter = max_iter, stop_early_exag_iter = stop_early_exag_iter, mom_switch_iter = mom_switch_iter)


def calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state, init = 'spectral', n_epochs = None, learning_rate  = 1.0, knn_indices = None, knn_dists = None):
	umap = UMAP(n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, spread = spread, random_state = random_state, \
		init = init, n_epochs = n_epochs, learning_rate = learning_rate, verbose = True)
	return umap.fit_transform(X, knn_indices = knn_indices, knn_dists = knn_dists)


def calc_force_directed_layout(W, file_name, n_jobs, target_change_per_node, target_steps, is3d, memory, random_state, init = None):
	input_graph_file = '{file_name}.net'.format(file_name = file_name)
	G = construct_graph(W)
	G.write(input_graph_file)
	print(input_graph_file + ' is written.')			

	output_coord_file = '{file_name}.coords.txt'.format(file_name = file_name)

	classpath = pkg_resources.resource_filename('scCloud', 'ext/forceatlas2.jar') + ':' + pkg_resources.resource_filename('scCloud', 'ext/gephi-toolkit-0.9.2-all.jar')
	command = ['java', '-Djava.awt.headless=true', '-Xmx{memory}g'.format(memory = memory), '-cp', classpath, \
				'kco.forceatlas2.Main', '--input', input_graph_file, '--output', file_name + '.coords', \
				'--nthreads', str(n_jobs), '--seed', str(random_state), '--targetChangePerNode', str(target_change_per_node), '--targetSteps', str(target_steps)]
	if not is3d:
		command.append('--2d')

	if init is not None:
		if not is3d:
			df = pd.DataFrame(data = init, columns = ['x', 'y'], index = range(1, init.shape[0] + 1))
		else:
			df = pd.DataFrame(data = init, columns = ['x', 'y', 'z'], index = range(1, init.shape[0] + 1))
		df.index.name = 'id'
		init_coord_file = '{file_name}.init.coords.txt'.format(file_name = file_name)
		df.to_csv(init_coord_file, sep = '\t', float_format = '%.2f')
		command.extend(['--coords', init_coord_file])
	print(command)
	check_call(command)
	print("Force-directed layout is generated.")

	fle_coords = pd.read_csv(output_coord_file, header = 0, index_col = 0, sep = '\t').values
	check_call(['rm', '-f', input_graph_file])
	check_call(['rm', '-f', output_coord_file])

	return fle_coords



def run_tsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, out_basis = 'tsne'):
	start = time.time()
	X_tsne = calc_tsne(data.obsm[rep_key], n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)
	data.obsm['X_' + out_basis] = X_tsne
	end = time.time()
	print("tSNE is calculated. Time spent = {:.2f}s.".format(end - start))


def run_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, out_basis = 'fitsne'):
	start = time.time()
	data.obsm['X_' + out_basis] = calc_fitsne(data.obsm[rep_key], n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state) 
	end = time.time()
	print("FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))


def run_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.5, spread = 1.0, random_state = 0, out_basis = 'umap'):
	start = time.time()
	assert 'knn_indices' in data.uns
	knn_indices = np.insert(data.uns['knn_indices'][:, 0 : n_neighbors - 1], 0, range(data.shape[0]), axis = 1)
	knn_dists = np.insert(data.uns['knn_distances'][:, 0 : n_neighbors - 1], 0, 0.0, axis = 1)
	data.obsm['X_' + out_basis] = calc_umap(data.obsm[rep_key], n_components, n_neighbors, min_dist, spread, random_state, knn_indices = knn_indices, knn_dists = knn_dists)
	end = time.time()
	print("UMAP is calculated. Time spent = {:.2f}s.".format(end - start))


def run_force_directed_layout(data, file_name, n_jobs, K = 50, target_change_per_node = 2.0, target_steps = 5000, is3d = False, memory = 8, random_state = 0, out_basis = 'fle'):
	start = time.time()
	realK = min(K - 1, data.uns['diffmap_knn_indices'].shape[1]) # K - 1: exclude self
	W = calculate_affinity_matrix(data.uns['diffmap_knn_indices'][:, 0:realK], data.uns['diffmap_knn_distances'][:, 0:realK])
	data.obsm['X_' + out_basis] = calc_force_directed_layout(W, file_name, n_jobs, target_change_per_node, target_steps, is3d, memory, random_state);
	end = time.time()
	print("Force-directed layout is calculated. Time spent = {:.2f}s.".format(end - start))



def run_net_tsne(data, rep_key, selected, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, net_alpha = 0.1, polish_learning_frac = 0.33, polish_n_iter = 150, out_basis = 'net_tsne'):
	start = time.time()
	X = data.obsm[rep_key][selected,:]
	X_tsne = calc_tsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)

	data.uns['X_' + out_basis + '_small'] = X_tsne
	data.obs['ds_selected'] = selected

	Y_init = np.zeros((data.shape[0], 2), dtype = np.float64)
	Y_init[selected,:] = X_tsne
	Y_init[~selected,:] = net_train_and_predict(X, X_tsne, data.obsm[rep_key][~selected,:], net_alpha, random_state, verbose = True)

	data.obsm['X_' + out_basis + '_pred'] = Y_init

	polish_learning_rate = polish_learning_frac * data.shape[0]
	data.obsm['X_' + out_basis] = calc_tsne(data.obsm[rep_key], n_jobs, n_components, perplexity, early_exaggeration, polish_learning_rate, random_state, \
		init = Y_init, n_iter = polish_n_iter, n_iter_early_exag = 0)
	end = time.time()
	print("Net tSNE is calculated. Time spent = {:.2f}s.".format(end - start))


def run_net_fitsne(data, rep_key, selected, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0, net_alpha = 0.1, polish_learning_frac = 0.5, polish_n_iter = 150, out_basis = 'net_fitsne'):
	start = time.time()
	X = data.obsm[rep_key][selected,:]
	X_fitsne = calc_fitsne(X, n_jobs, n_components, perplexity, early_exaggeration, learning_rate, random_state)

	data.uns['X_' + out_basis + '_small'] = X_fitsne
	data.obs['ds_selected'] = selected

	Y_init = np.zeros((data.shape[0], 2), dtype = np.float64)
	Y_init[selected,:] = X_fitsne
	Y_init[~selected,:] = net_train_and_predict(X, X_fitsne, data.obsm[rep_key][~selected,:], net_alpha, random_state, verbose = True)

	data.obsm['X_' + out_basis + '_pred'] = Y_init

	polish_learning_rate = polish_learning_frac * data.shape[0]
	data.obsm['X_' + out_basis] = calc_fitsne(data.obsm[rep_key], n_jobs, n_components, perplexity, early_exaggeration, polish_learning_rate, random_state, 
		initialization = Y_init, max_iter = polish_n_iter, stop_early_exag_iter = 0, mom_switch_iter = 0)
	end = time.time()
	print("Net FItSNE is calculated. Time spent = {:.2f}s.".format(end - start))


def run_net_umap(data, rep_key, selected, n_jobs, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 0, ds_full_speed = False, net_alpha = 0.1, polish_n_epochs = 30, polish_learning_rate = 10.0, out_basis = 'net_umap'):
	start = time.time()
	X = data.obsm[rep_key][selected,:]

	if 'ds_knn_indices' not in data.uns:
		indices, distances = calculate_nearest_neighbors(X, n_jobs, K = n_neighbors, random_state = random_state, full_speed = ds_full_speed)
		data.uns['ds_knn_indices'] = indices
		data.uns['ds_knn_distances'] = distances
	knn_indices = np.insert(data.uns['ds_knn_indices'][:, 0 : n_neighbors - 1], 0, range(X.shape[0]), axis = 1)
	knn_dists = np.insert(data.uns['ds_knn_distances'][:, 0 : n_neighbors - 1], 0, 0.0, axis = 1)

	X_umap = calc_umap(X, n_components, n_neighbors, min_dist, spread, random_state, knn_indices = knn_indices, knn_dists = knn_dists)

	data.uns['X_' + out_basis + '_small'] = X_umap
	data.obs['ds_selected'] = selected
	
	Y_init = np.zeros((data.shape[0], 2), dtype = np.float64)
	Y_init[selected,:] = X_umap
	Y_init[~selected,:] = net_train_and_predict(X, X_umap, data.obsm[rep_key][~selected,:], net_alpha, random_state, verbose = True)

	data.obsm['X_' + out_basis + '_pred'] = Y_init

	assert 'knn_indices' in data.uns
	knn_indices = np.insert(data.uns['knn_indices'][:, 0 : n_neighbors - 1], 0, range(data.shape[0]), axis = 1)
	knn_dists = np.insert(data.uns['knn_distances'][:, 0 : n_neighbors - 1], 0, 0.0, axis = 1)

	data.obsm['X_' + out_basis] = calc_umap(data.obsm[rep_key], n_components, n_neighbors, min_dist, spread, random_state, \
		init = Y_init, n_epochs = polish_n_epochs, learning_rate = polish_learning_rate, knn_indices = knn_indices, knn_dists = knn_dists)
	end = time.time()
	print("Net UMAP is calculated. Time spent = {:.2f}s.".format(end - start))

def run_net_fle(data, selected, file_name, n_jobs, K = 50, target_change_per_node = 2.0, target_steps = 10000, is3d = False, memory = 8, random_state = 0, \
	ds_full_speed = False, net_alpha = 0.1, polish_target_steps = 1500, out_basis = 'net_fle'):
	start = time.time()
	X = data.obsm['X_diffmap'][selected,:]
	
	if 'diffmap_ds_knn_indices' in data.uns:
		indices = data.uns['diffmap_ds_knn_indices']
		distances = data.uns['diffmap_ds_knn_distances']
	else:
		indices, distances = calculate_nearest_neighbors(X, n_jobs, K = K, random_state = random_state, full_speed = ds_full_speed)
		data.uns['diffmap_ds_knn_indices'] = indices
		data.uns['diffmap_ds_knn_distances'] = distances

	W = calculate_affinity_matrix(indices, distances)
	X_fle = calc_force_directed_layout(W, file_name + '.small', n_jobs, target_change_per_node, target_steps, is3d, memory, random_state);

	data.uns['X_' + out_basis + '_small'] = X_fle
	data.obs['ds_diffmap_selected'] = selected

	Y_init = np.zeros((data.shape[0], 2), dtype = np.float64)
	Y_init[selected,:] = X_fle
	Y_init[~selected,:] = net_train_and_predict(X, X_fle, data.obsm['X_diffmap'][~selected,:], net_alpha, random_state, verbose = True)

	data.obsm['X_' + out_basis + '_pred'] = Y_init

	realK = min(K - 1, data.uns['diffmap_knn_indices'].shape[1]) # K - 1: exclude self
	W = calculate_affinity_matrix(data.uns['diffmap_knn_indices'][:, 0:realK], data.uns['diffmap_knn_distances'][:, 0:realK])
	data.obsm['X_' + out_basis] = calc_force_directed_layout(W, file_name, n_jobs, target_change_per_node, polish_target_steps, is3d, memory, random_state, init = Y_init);
	end = time.time()

	print("Net FLE is calculated. Time spent = {:.2f}s.".format(end - start))
