import time
import numpy as np
import pandas as pd

from scipy.sparse import issparse
from scipy.stats import entropy, chi2
from sklearn.neighbors import NearestNeighbors



def calculate_nearest_neighbors(X, num_threads, method = 'hnsw', K = 100, M = 20, efC = 200, efS = 200, random_state = 0, full_speed = False):
	"""X is the sample by feature matrix, could be either dense or sparse"""

	start_time = time.time()

	nsample = X.shape[0]
	
	if nsample <= 1000:
		method = 'sklearn'

	if nsample < K:
		print("Warning: in calculate_nearest_neighbors, number of samples = {} < K = {}!\n Set K to {}.".format(nsample, K, nsample))
		K = nsample

	knn_index = None
	if method == 'hnsw':
		import hnswlib
		assert not issparse(X)
		# Build hnsw index
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



def select_cells(distances, frac, K = 25, alpha = 1.0, random_state = 0):
	start_time = time.time()

	nsample = distances.shape[0]

	if K > distances.shape[1]:
		print("Warning: in select_cells, K = {} > the number of calculated nearest neighbors!\nSet K to {}".format(K, distances.shape[1]))
		K = distances.shape[1]

	probs = np.zeros(nsample)
	if alpha == 0.0:
		probs[:] = 1.0 # uniform
	elif alpha == 1.0:
		probs[:] = distances[:, K - 1]
	else:
		probs[:] = distances[:, K - 1] ** alpha
	probs /= probs.sum()

	np.random.seed(random_state)
	selected = np.zeros(nsample, dtype = bool)
	selected[np.random.choice(nsample, size = int(nsample * frac), replace = False, p = probs)] = True
	
	end_time = time.time()
	print('select_cells finished. Time spent = {:.2}s.'.format(end_time - start_time))
	
	return selected



def calc_kBET_for_one_datapoint(pos, attr_values, knn_indices, ideal_dist, K):
	indices = np.append(knn_indices[pos], [pos])
	df = ideal_dist.size - 1

	observed_counts = pd.Series(attr_values[indices]).value_counts(sort = False).values
	expected_counts = ideal_dist * K
	static = np.sum(np.divide(np.square(np.subtract(observed_counts, expected_counts)), expected_counts))
	p_value = 1 - chi2.cdf(static, df)
	return (static, p_value)


def calc_kBET(data, attr, knn_keyword = 'knn', K = 25, n_jobs = 1, temp_folder = None):
	"""
	This kBET metric is based on paper "A test metric for assessing single-cell RNA-seq batch correction" [M. Büttner, et al.] in Nature Methods, 2018.

	:return:
		static_mean: average chi-square static over all the data points.
		pvalue_mean: average p-value over all the data points.
	"""
	assert attr in data.obs and data.obs[attr].dtype.name == 'category'

	from joblib import Parallel, delayed

	ideal_dist = data.obs[attr].value_counts(normalize = True, sort = False).values # ideal no batch effect distribution
	nsample = data.shape[0]	
	nbatch = ideal_dist.size

	attr_values = data.obs[attr].values.copy()
	attr_values.categories = range(nbatch)

	knn_indices = data.uns[knn_keyword + '_indices'][:, 0 : K - 1]
	kBET_arr = np.array(Parallel(n_jobs = n_jobs, max_nbytes = 1e7, temp_folder = temp_folder)(delayed(calc_kBET_for_one_datapoint)(i, attr_values, knn_indices, ideal_dist, K) for i in range(nsample)))

	res = kBET_arr.mean(axis = 0)
	static_mean = res[0]
	pvalue_mean = res[1]

	return (static_mean, pvalue_mean)




# def calc_JSD(P, Q):
# 	M = (P + Q) / 2
# 	return (entropy(P, M, base = 2) + entropy(Q, M, base = 2)) / 2.0

# def calc_kBJSD_for_one_datapoint(pos, attr_values, knn_indices, ideal_dist):
# 	idx = np.append(knn_indices[pos], [pos])
# 	empirical_dist = pd.Series(attr_values[idx]).value_counts(normalize = True, sort = False).values
# 	return calc_JSD(ideal_dist, empirical_dist)

# def calc_kBJSD(data, attr, knn_keyword = 'knn', K = 25, n_jobs = 1, temp_folder = None):
# 	assert attr in data.obs and data.obs[attr].dtype.name == 'category'

# 	from joblib import Parallel, delayed

# 	ideal_dist = data.obs[attr].value_counts(normalize = True, sort = False).values # ideal no batch effect distribution
# 	nsample = data.shape[0]	
# 	nbatch = ideal_dist.size

# 	attr_values = data.obs[attr].values.copy()
# 	attr_values.categories = range(nbatch)
# 	#attr_values = attr_values.astype(int)

# 	knn_indices = data.uns[knn_keyword + '_indices'][:, 0 : K - 1]
# 	kBJSD_arr = np.array(Parallel(n_jobs = n_jobs, max_nbytes = 1e7, temp_folder = temp_folder)(delayed(calc_kBJSD_for_one_datapoint)(i, attr_values, knn_indices, ideal_dist) for i in range(nsample)))

# 	return kBJSD_arr.mean()
