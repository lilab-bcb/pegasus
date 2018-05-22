import time
import numpy as np
import pandas as pd

import igraph
import louvain
import hdbscan
from sklearn.cluster import KMeans
from natsort import natsorted


def construct_graph(W, directed = True):
	s, t = W.nonzero()
	w = W[s, t].A1
	G = igraph.Graph(directed = directed)
	G.add_vertices(W.shape[0])
	G.add_edges(zip(s, t))
	G.es['weight'] = w
	assert G.vcount() == W.shape[0]

	return G

def run_louvain(data, resolution = 1.3, random_state = 0):
	start = time.time()
	louvain.set_rng_seed(random_state)
	G = construct_graph(data.uns['W_norm'])
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
	categories = natsorted(label_map.values())
	data.obs['hdbscan_labels'] = pd.Categorical(values = labels, categories = categories)

	soft_clusters = np.argmax(hdbscan.all_points_membership_vectors(clusterer), axis = 1)
	labels[noise_idx] = f_trans(soft_clusters[noise_idx])
	data.obs['hdbscan_labels_soft'] = pd.Categorical(values = labels, categories = categories)

	end = time.time()
	print("HDBSCAN clustering is done. Time spent = {:.2f}s.".format(end - start))



def run_kmeans(data, rep_key, n_clusters, n_jobs = 1):
	start = time.time()
	km = KMeans(n_clusters = 20, n_jobs = 10)
	km.fit(data.obsm[rep_key].astype('float64'))
	ids, counts = np.unique(km.labels_, return_counts = True)
	label_map = dict(zip(ids[np.argsort(counts)[::-1]], [str(x + 1) for x in range(len(counts))]))
	f_trans = np.vectorize(lambda x: label_map[x])
	labels = f_trans(km.labels_)
	categories = natsorted(label_map.values())
	data.obs['kmeans_labels'] = pd.Categorical(values = labels, categories = categories)

	end = time.time()
	print("Time spent = {:.2f}".format(end - start))

