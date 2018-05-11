import time
import numpy as np
import pandas as pd

from MulticoreTSNE import MulticoreTSNE as TSNE
from umap import UMAP
from fitsne import FItSNE

def run_tsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, learning_rate = 1000, random_state = 0):
	start = time.time()
	tsne = TSNE(n_jobs = n_jobs, n_components = n_components, perplexity = perplexity, early_exaggeration = early_exaggeration, learning_rate = learning_rate, random_state = random_state)
	X = data.obsm[rep_key].astype('float64')
	X_tsne = tsne.fit_transform(X)
	data.obsm['X_tsne'] = X_tsne
	end = time.time()
	print("TSNE is done. Time spent = {:.2f}s.".format(end - start))

def run_fitsne(data, rep_key, n_jobs, n_components = 2, perplexity = 30, early_exaggeration = 12, random_state = 0):
	start = time.time()
	X = data.obsm[rep_key].astype('float64').copy(order = 'C')
	X_fitsne = FItSNE(X, nthreads = n_jobs, no_dims = n_components, perplexity = perplexity, rand_seed = (random_state if random_state is not None else -1), early_exag_coeff = early_exaggeration)
	data.obsm['X_fitsne'] = X_fitsne
	end = time.time()
	print("FItSNE is done. Time spent = {:.2f}s.".format(end - start))

def run_umap(data, rep_key, n_components = 2, n_neighbors = 15, min_dist = 0.1, spread = 1.0, random_state = 0):
	start = time.time()
	umap = UMAP(n_components = n_components, n_neighbors = n_neighbors, min_dist = min_dist, spread = spread, random_state = random_state)
	X = data.obsm[rep_key].astype('float64')
	X_umap = umap.fit_transform(X)
	data.obsm['X_umap'] = X_umap	
	end = time.time()
	print("UMAP is done. Time spent = {:.2f}s.".format(end - start))

def run_force_directed_layout(data, rep_key, n_jobs):
	start = time.time()
	# important codes
	end = time.time()
	print("Force-directed layout is done. Time spent = {:.2f}s.".format(end - start))
