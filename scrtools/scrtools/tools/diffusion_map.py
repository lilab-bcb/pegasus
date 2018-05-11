import time 
import numpy as np
import pandas as pd
import nmslib 

from scipy.sparse import issparse, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from sklearn.utils.extmath import randomized_svd
from sklearn.neighbors import NearestNeighbors

def get_symmetric_matrix(csr_mat):
	tp_mat = csr_mat.transpose().tocsr()
	sym_mat = csr_mat + tp_mat
	sym_mat.sort_indices()

	idx_mat = (csr_mat != 0).astype(int) + (tp_mat != 0).astype(int)
	idx_mat.sort_indices()
	idx = idx_mat.data == 2

	sym_mat.data[idx] /= 2.0
	return sym_mat

def calculate_affinity_matrix(X, num_threads, method = 'nmslib', K = 100, M = 15, efC = 100, efS = 100):
	"""X is the sample by feature matrix, could be either dense or sparse"""
	nsample = X.shape[0]
	
	if method == 'nmslib':
		if issparse(X):
			knn_index = nmslib.init(method = 'hnsw', space = 'l2_sparse', data_type = nmslib.DataType.SPARSE_VECTOR)
		else:
			knn_index = nmslib.init(method = 'hnsw', space = 'l2', data_type = nmslib.DataType.DENSE_VECTOR)

		knn_index.addDataPointBatch(X)
		start = time.time()
		knn_index.createIndex({'M': M, 'indexThreadQty': num_threads, 'efConstruction': efC}) 
		end = time.time() 
		print("Indexing time = {0:.6f}".format(end - start))
		
		knn_index.setQueryTimeParams({'efSearch': efS})
		start = time.time() 
		results = knn_index.knnQueryBatch(X, k = K, num_threads = num_threads)
		end = time.time() 
		print("kNN time total={0:.6f} (sec), per query={1:.6f} (sec), per query adjusted for thread number={2:.6f} (sec)".format(end - start, float(end - start) / nsample, num_threads * (float(end - start) / nsample)))

		K = K - 1 # eliminate the first neighbor, which is the node itself
		indices = np.zeros((nsample, K), dtype = int)
		distances = np.zeros((nsample, K))
		for i in range(nsample):
			if results[i][0][0] == i:
				indices[i,:] = results[i][0][1:]
				distances[i,:] = results[i][1][1:]
			else:
				indices[i,:] = results[i][0][:K]
				distances[i,:] = results[i][1][:K]
	else:
		assert method == 'sklearn'
		K = K - 1
		knn = NearestNeighbors(n_neighbors = K, n_jobs = num_threads)
		knn.fit(X)
		distances, indices = knn.kneighbors()


	# calculate sigma, important to use median here!
	sigmas = np.median(distances, axis = 1)
	sigmas_sq = np.square(sigmas)

	# calculate local-scaled kernel
	for i in range(nsample):
		numers = 2.0 * sigmas[i] * sigmas[indices[i,:]]
		denoms = sigmas_sq[i] + sigmas_sq[indices[i,:]]
		distances[i,:] = np.sqrt(numers / denoms) * np.exp(-np.square(distances[i,:]) / denoms)
		
	W = csr_matrix((distances.ravel(), (np.repeat(range(nsample), K), indices.ravel())))
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

def calculate_diffusion_map(W, n_dc = 100, alpha = 0.5, solver = 'randomized', random_state = None):
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

	U_pt = U * S_new #symmetric pseudo component
	Phi_pt = Phi * S_new #asym pseudo component

	end = time.time()
	print("Calculating diffusion map is done.")

	return Phi_pt, U_pt, S, W_norm

def run_diffmap(data, rep_key, n_jobs = 1, n_components = 100, alpha = 0.5, K = 100, random_state = 0, knn_method = 'nmslib', eigen_solver = 'randomized', M = 15, efC = 100, efS = 100):
	start = time.time()
	W = calculate_affinity_matrix(data.obsm[rep_key], n_jobs, method = knn_method, K = K, M = M, efC = efC, efS = efS)
	Phi_pt, U_pt, S, W_norm = calculate_diffusion_map(W, n_dc = n_components, alpha = alpha, solver = eigen_solver, random_state = random_state)
	data.uns['W'] = W
	data.uns['W_norm'] = W_norm
	data.uns['diffmap_evals'] = S
	data.obsm['X_diffmap'] = Phi_pt
	data.obsm['X_diffmap_sym'] = U_pt
	end = time.time()
	print("run_diffmap finished. Time spent = {:.2f}s.".format(end - start))
