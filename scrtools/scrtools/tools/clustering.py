import time
import numpy as np
import pandas as pd

import igraph
import louvain
import hdbscan
from sklearn.cluster import KMeans
from natsort import natsorted
import threading

import ctypes
import ctypes.util

def construct_graph(W, directed = True):
	s, t = W.nonzero()
	w = W[s, t].A1
	G = igraph.Graph(directed = directed)
	G.add_vertices(W.shape[0])
	G.add_edges(zip(s, t))
	G.es['weight'] = w
	assert G.vcount() == W.shape[0]

	return G

def run_louvain(data, affinity = 'W_norm', resolution = 1.3, random_state = 0):
	start = time.time()
	louvain.set_rng_seed(random_state)
	G = construct_graph(data.uns[affinity])
	partition = louvain.find_partition(G, louvain.RBConfigurationVertexPartition, resolution_parameter = resolution)
	labels = np.array([str(x + 1) for x in partition.membership])
	categories = natsorted(np.unique(labels))
	data.obs['louvain_labels'] = pd.Categorical(values = labels, categories = categories)
	end = time.time()
	print("Louvain clustering is done. Time spent = {:.2f}s.".format(end - start))



def run_hdbscan(data, rep_key, n_jobs = 1, min_cluster_size = 50, min_samples = 25):
	start = time.time()
	clusterer = hdbscan.HDBSCAN(core_dist_n_jobs = n_jobs, min_cluster_size = min_cluster_size, min_samples = min_samples, prediction_data = True)
	clusterer.fit(data.obsm[rep_key].astype('float64'))
	
	noise_idx = clusterer.labels_ < 0
	ids, counts = np.unique(clusterer.labels_[~noise_idx], return_counts = True)
	label_map = dict(zip(ids[np.argsort(counts)[::-1]], [str(x + 1) for x in range(len(counts))]))
	f_trans = np.vectorize(lambda x: label_map.get(x, 'noise'))
	
	labels = f_trans(clusterer.labels_)
	categories = natsorted(list(label_map.values()) + ['noise'])
	data.obs['hdbscan_labels'] = pd.Categorical(values = labels, categories = categories)

	soft_clusters = np.argmax(hdbscan.all_points_membership_vectors(clusterer), axis = 1)
	labels[noise_idx] = f_trans(soft_clusters[noise_idx])
	data.obs['hdbscan_labels_soft'] = pd.Categorical(values = labels, categories = categories)

	end = time.time()
	print("HDBSCAN clustering is done. Time spent = {:.2f}s.".format(end - start))



def set_numpy_thread(number):
	old_n = 0
	openblas_loc = ctypes.util.find_library('openblas')
	if openblas_loc is not None:
		openblas_lib = ctypes.cdll.LoadLibrary(openblas_loc)
		old_n = openblas_lib.openblas_get_num_threads()
		openblas_lib.openblas_set_num_threads(number)
	else:
		mkl_loc = ctypes.util.find_library('mkl_rt')
		if mkl_loc is not None:
			mkl_lib = ctypes.cdll.LoadLibrary(mkl_loc)
			old_n = mkl_lib.mkl_get_max_threads()
			mkl_lib.mkl_set_num_threads(ctypes.byref(ctypes.c_int(number)))
	return old_n


def run_kmeans(data, rep_key, n_clusters, n_init = 10, n_jobs = 1, random_state = 0):
	start = time.time()

	old_n = set_numpy_thread(1)
	km = KMeans(n_clusters = n_clusters, n_init = n_init, n_jobs = n_jobs, random_state = random_state)
	km.fit(data.obsm[rep_key].astype('float64'))
	set_numpy_thread(old_n)

	ids, counts = np.unique(km.labels_, return_counts = True)
	label_map = dict(zip(ids[np.argsort(counts)[::-1]], [str(x + 1) for x in range(len(counts))]))
	f_trans = np.vectorize(lambda x: label_map[x])
	
	labels = f_trans(km.labels_)
	categories = natsorted(label_map.values())
	data.obs['kmeans_labels'] = pd.Categorical(values = labels, categories = categories)

	end = time.time()
	print("Spectral clustering is done. Time spent = {:.2f}s".format(end - start))


def run_one_instance_of_kmeans(thread_no, results, n_init, n_clusters, n_jobs, X, seeds):
	results[thread_no] = []
	for i in range(thread_no, n_init, n_jobs):
		km = KMeans(n_clusters = n_clusters, n_init = 1, n_jobs = 1, random_state = seeds[i])
		km.fit(X)
		results[thread_no].append(km.labels_)


def run_approximated_louvain(data, rep_key, n_jobs = 1, resolution = 1.3, random_state = 0, n_clusters = 30, n_init = 20):
	start = time.time()

	X = data.obsm[rep_key].astype('float64')
	np.random.seed(random_state)
	seeds = np.random.randint(np.iinfo(np.int32).max, size = n_init)
	
	old_n = set_numpy_thread(1)

	threads = [None] * n_jobs
	results = [None] * n_jobs

	for i in range(n_jobs):
		t = threading.Thread(target=run_one_instance_of_kmeans, args=(i, results, n_init, n_clusters, n_jobs, X, seeds))
		threads[i] = t
		t.start()

	for i in range(n_jobs):
		threads[i].join()

	set_numpy_thread(old_n)
	
	labels = list(zip(*[x for y in results for x in y]))
	uniqs = np.unique(labels, axis = 0)
	transfer_dict = {tuple(k):v for k, v in zip(uniqs, range(uniqs.shape[0]))}
	labels = [transfer_dict[x] for x in labels]

	louvain.set_rng_seed(random_state)
	G = construct_graph(data.uns['W_norm'])
	partition = louvain.RBConfigurationVertexPartition(G, resolution_parameter = resolution, initial_membership = labels)

	partition_agg = partition.aggregate_partition()
	optimiser = louvain.Optimiser()
	diff = optimiser.optimise_partition(partition_agg)
	partition.from_coarse_partition(partition_agg)

	labels = np.array([str(x + 1) for x in partition.membership])
	categories = natsorted(np.unique(labels))
	data.obs['approx_louvain_labels'] = pd.Categorical(values = labels, categories = categories)

	end = time.time()
	print("Approximated Louvain clustering is done. Time spent = {:.2f}s.".format(end - start))
