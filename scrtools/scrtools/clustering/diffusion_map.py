import numpy as np
import pandas as pd
import nmslib 
import time 

from scipy.sparse import issparse, csr_matrix, diags
from scipy.sparse.csgraph import connected_components
from sklearn.utils.extmath import randomized_svd

def calculate_symmetric_knn_graph(X, num_threads, K = 100, M = 15, efC = 100, efS = 100, gamma = None):
	"""X is the sample by feature matrix, could be either dense or sparse; gamma default is 1 / n_features"""
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
	nsample = X.shape[0]
	start = time.time() 
	results = knn_index.knnQueryBatch(X, k = K, num_threads = num_threads)
	end = time.time() 
	print("kNN time total={0:.6f} (sec), per query={1:.6f} (sec), per query adjusted for thread number={2:.6f} (sec)".format(end - start, float(end - start) / nsample, num_threads * (float(end - start) / nsample)))

	# Construct sparse symmetric kNN with Gaussian heat kernel
	ij = np.zeros((2, nsample * K), dtype = int)
	ij[0,] = np.repeat(range(nsample), K)
	ij[1,] = np.concatenate([x[0] for x in results])
	data = np.concatenate([x[1] for x in results])

	for i in range(nsample):
		if ij[1, i * K] != i:
			ij[1, i * K + K - 1] = i
			data[i * K + K - 1] = 0.0

	if gamma is None:
		gamma = 1.0 / X.shape[1]
	data = np.exp(-gamma * (data ** 2))

	csr_mat = csr_matrix((data, ij))
	tp_mat = csr_mat.transpose().tocsr()
	all_mat = csr_mat + tp_mat
	all_mat.sort_indices()

	idx_mat = (csr_mat != 0).astype(int) + (tp_mat != 0).astype(int)
	idx_mat.sort_indices()
	idx = idx_mat.data == 2

	all_mat.data[idx] /= 2.0

	return all_mat

def calculate_diffusion_map(W, n_components = 15, sigma = 0.01, random_state = None, return_eigen_vectors = True):
	assert issparse(W)

	start = time.time()
	nc, labels = connected_components(W, directed = True, connection = 'strong')
	end = time.time()
	print("Finding connected components time = {0:.6f}".format(end - start))

	assert nc == 1

	diag = W.sum(axis = 1).A1
	diag_half = np.sqrt(diag)
	W_norm = W.tocoo(copy = True)
	W_norm.data /= diag_half[W_norm.row]
	W_norm.data /= diag_half[W_norm.col]

	start = time.time()
	U, S, V = randomized_svd(W_norm, n_components = n_components, random_state = random_state)
	end = time.time()
	print("Randomized SVD time = {0:.6f}".format(end - start))

	Phi = U / diag_half[:, None]

	t = int(np.log(sigma) / np.log(S[-1])) # diffusion time
	Phi_t = Phi * (S ** t)[None, :]

	results = [Phi_t, t]
	if return_eigen_vectors:
		results.append(Phi)

	return results
