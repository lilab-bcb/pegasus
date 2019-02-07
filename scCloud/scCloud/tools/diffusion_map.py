import time 
import numpy as np
import pandas as pd

from scipy.sparse import issparse, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from sklearn.utils.extmath import randomized_svd
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances

from . import calculate_nearest_neighbors



def get_symmetric_matrix(csr_mat):
	tp_mat = csr_mat.transpose().tocsr()
	sym_mat = csr_mat + tp_mat
	sym_mat.sort_indices()

	idx_mat = (csr_mat != 0).astype(int) + (tp_mat != 0).astype(int)
	idx_mat.sort_indices()
	idx = idx_mat.data == 2

	sym_mat.data[idx] /= 2.0
	return sym_mat

def calculate_affinity_matrix(indices, distances):
	nsample = indices.shape[0]
	K = indices.shape[1]
	# calculate sigma, important to use median here!
	sigmas = np.median(distances, axis = 1)
	sigmas_sq = np.square(sigmas)

	# calculate local-scaled kernel
	for i in range(nsample):
		numers = 2.0 * sigmas[i] * sigmas[indices[i,:]]
		denoms = sigmas_sq[i] + sigmas_sq[indices[i,:]]
		distances[i,:] = np.sqrt(numers / denoms) * np.exp(-np.square(distances[i,:]) / denoms)

	W = csr_matrix((distances.ravel(), (np.repeat(range(nsample), K), indices.ravel())), shape = (nsample, nsample))
	W = get_symmetric_matrix(W)

	# density normalization
	z = W.sum(axis = 1).A1
	W = W.tocoo()
	W.data /= z[W.row]
	W.data /= z[W.col]
	W = W.tocsr()

	print("Constructing affinity matrix is done.")

	return W

def calculate_normalized_affinity(W):
	diag = W.sum(axis = 1).A1
	diag_half = np.sqrt(diag)
	W_norm = W.tocoo(copy = True)
	W_norm.data /= diag_half[W_norm.row]
	W_norm.data /= diag_half[W_norm.col]
	W_norm = W_norm.tocsr()

	return W_norm, diag, diag_half

def calculate_diffusion_map(W, n_dc = 100, alpha = 0.5, solver = 'randomized', random_state = 0):
	assert issparse(W)

	start = time.time()
	nc, labels = connected_components(W, directed = True, connection = 'strong')
	print("Calculating connected components is done.")

	assert nc == 1

	W_norm, diag, diag_half = calculate_normalized_affinity(W)
	print("Calculating normalized affinity matrix is done.")

	if solver == 'randomized':
		U, S, VT = randomized_svd(W_norm, n_components = n_dc, random_state = random_state)
	else:
		assert solver == 'eigsh'
		S, U = eigsh(W_norm, k = n_dc)
		S = S[::-1]
		U = U[:, ::-1]

	# remove the first eigen value and vector
	S = S[1:]
	U = U[:, 1:]

	Phi = U / diag_half[:, np.newaxis]
	S_new = S / (1 - alpha * S)

	U_df = U * S #symmetric diffusion component
	Phi_pt = Phi * S_new #asym pseudo component

	end = time.time()
	print("Calculating diffusion map is done.")

	return Phi_pt, U_df, S, W_norm

def reduce_diffmap_to_3d(Phi_pt, random_state = 0):
	pca = PCA(n_components = 3, random_state = random_state)
	Phi_reduced = pca.fit_transform(Phi_pt)
	print("Reduce diffmap to 3D is done.")

	return Phi_reduced

def run_diffmap(data, rep_key, n_jobs = 1, n_components = 100, alpha = 0.5, K = 100, random_state = 0, knn_method = 'hnsw', eigen_solver = 'randomized', M = 20, efC = 200, efS = 200, full_speed = False):
	start = time.time()

	indices, distances, knn_index = calculate_nearest_neighbors(data.obsm[rep_key], n_jobs, method = knn_method, \
		K = K, M = M, efC = efC, efS = efS, random_state = random_state, full_speed = full_speed)
	if knn_index is not None:
		data.uns['knn'] = knn_index
		data.uns['knn_dim'] = data.obsm[rep_key].shape[1]
	W = calculate_affinity_matrix(indices, distances)

	Phi_pt, U_df, S, W_norm = calculate_diffusion_map(W, n_dc = n_components, alpha = alpha, solver = eigen_solver, random_state = random_state)

	Phi_reduced = reduce_diffmap_to_3d(Phi_pt, random_state = random_state)

	indices, distances, knn_index = calculate_nearest_neighbors(Phi_pt, n_jobs, method = knn_method, \
		K = K, M = M, efC = efC, efS = efS, random_state = random_state, full_speed = full_speed)
	if knn_index is not None:
		data.uns['diffmap_knn'] = knn_index
		data.uns['diffmap_knn_dim'] = Phi_pt.shape[1]
	W_diffmap = calculate_affinity_matrix(indices, distances)
	
	# W_diffmap_norm, diag_tmp, diag_half_tmp = calculate_normalized_affinity(W_diffmap)

	data.uns['W'] = W
	data.uns['W_norm'] = W_norm
	data.uns['diffmap_evals'] = S
	data.uns['W_diffmap'] = W_diffmap
	# data.uns['W_diffmap_norm'] = W_diffmap_norm

	data.obsm['X_diffmap'] = Phi_pt
	data.obsm['X_dmnorm'] = U_df
	data.obsm['X_diffmap_pca'] = Phi_reduced
	end = time.time()
	print("run_diffmap finished. Time spent = {:.2f}s.".format(end - start))



def run_pseudotime_calculation(data, roots):
	start = time.time()
	data.uns['roots'] = roots
	mask = np.isin(data.obs_names, data.uns['roots'])
	distances = np.mean(euclidean_distances(data.obsm['X_diffmap'][mask, :], data.obsm['X_diffmap']), axis = 0)
	dmin = distances.min()
	dmax = distances.max()
	data.obs['pseudotime'] = (distances - dmin) / (dmax - dmin)
	end = time.time()
	print("run_pseudotime_calculation finished. Time spent = {:.2f}s".format(end - start))
