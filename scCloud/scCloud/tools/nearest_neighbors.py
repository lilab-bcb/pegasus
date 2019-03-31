import time
import numpy as np


from scipy.sparse import issparse
from sklearn.neighbors import NearestNeighbors



def calculate_nearest_neighbors(X, num_threads, method = 'hnsw', hnsw_index = None, K = 100, M = 20, efC = 200, efS = 200, random_state = 0, full_speed = False):
	"""X is the sample by feature matrix, could be either dense or sparse"""

	start_time = time.time()

	nsample = X.shape[0]
	if nsample < 500:
		method = 'sklearn'

	knn_index = None

	if method == 'hnsw':
		import hnswlib
		assert not issparse(X)
		knn_index = hnsw_index

		# Build hnsw index
		if knn_index is None:
			knn_index = hnswlib.Index(space = 'l2', dim = X.shape[1])
			knn_index.init_index(max_elements = nsample, ef_construction = efC, M = M, random_seed = random_state)
			knn_index.set_num_threads(num_threads if full_speed else 1)
			knn_index.add_items(X)

		# KNN query
		knn_index.set_ef(efS)
		knn_index.set_num_threads(num_threads)
		indices, distances = knn_index.knn_query(X, k = K)
		# eliminate the first neighbor, which is the node itself
		for i in range(nsample):
			if indices[i, 0] != i:
				indices[i, 1:] = indices[i, 0:-1]
				distances[i, 1:] = distances[i, 0:-1]
		indices = indices[:, 1:].astype(int)
		distances = np.sqrt(distances[:, 1:])
	else:
		assert method == 'sklearn'
		knn = NearestNeighbors(n_neighbors = K - 1, n_jobs = num_threads) # eliminate the first neighbor, which is the node itself
		knn.fit(X)
		distances, indices = knn.kneighbors()

	end_time = time.time()
	print("Nearest neighbor search is finished in {:.2f}s.".format(end_time - start_time))

	return (indices, distances)



def select_cells(indices, first_K, random_state = 0):
	nsample = indices.shape[0]
	K = indices.shape[1]
	np.random.seed(random_state)
	orders = np.random.permutation(nsample)
	selected = np.zeros(nsample, dtype = bool)
	covered = np.zeros(nsample, dtype = bool)
	for cell_id in orders:
		if not covered[cell_id]:
			covered[cell_id] = True
			selected[cell_id] = True
			for i in range(first_K):
				covered[indices[cell_id, i]] = True
	return selected



