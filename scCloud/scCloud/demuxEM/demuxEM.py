import numpy as np
import pandas as pd
import time
from natsort import natsorted

import multiprocessing

# global variables
pvec = None
alpha = 0.0
alpha_noise = 1.0
tol = 1e-6

def estimate_probs(arr, pvec, alpha, alpha_noise, tol):
	probs = np.zeros(pvec.size + 1)
	old_probs = np.zeros(pvec.size + 1)
	z = np.zeros(pvec.size + 1)
	noise = pvec.size
	# Estimate MLE without Generalized Dirichlet prior
	probs_mle = arr / arr.sum()
	probs[noise] = (probs_mle / pvec).min() + 0.01
	probs[:-1] = np.maximum(probs_mle - probs[noise] * pvec, 0.01)
	probs = probs / probs.sum()

	# EM algorithm
	i = 0
	eps = 1.0
	while eps > tol:
		i += 1
		old_probs[:] = probs[:]
		# E step
		z[:-1] = alpha - 1.0
		z[noise] = alpha_noise - 1.0
		for j in range(pvec.size):
			if arr[j] > 0:
				p = probs[j] / (probs[noise] * pvec[j] + probs[j])
				z[j] += arr[j] * p
				z[noise] += arr[j] * (1.0 - p)
		# M step
		idx = z > 0.0
		probs[idx] = z[idx] / z[idx].sum()
		probs[~idx] = 0.0
		eps = np.linalg.norm(probs - old_probs, ord = 1)
		# print ("i = {}, eps = {:.2g}.".format(i, eps))

	return probs

def estimate_probs_wrapper(arr):
	return estimate_probs(arr, pvec, alpha, alpha_noise, tol)

def get_droplet_info(arr):
	res = [str(x + 1) for x in arr.nonzero()[0]]
	return ('singlet' if len(res) == 1 else 'doublet', ','.join(res))

def demultiplex(data, adt, unknown = 0.8, alpha_value = 0.0, n_threads = 1):
	global pvec, alpha

	start = time.time()

	alpha = alpha_value # set alpha parameter for real samples
	
	idx_df = data.obs_names.isin(adt.obs_names)
	if idx_df.sum() < data.shape[0]:
		nzero = data.shape[0] - idx_df.sum()
		print("Warning: {} cells do not have ADTs, percentage = {:.2f}%.".format(nzero, nzero * 100.0 / data.shape[0]))

	adt_small = adt[data.obs_names[idx_df],].X.toarray()
	value = np.median(adt_small.sum(axis = 1))
	tmp = adt[adt.obs_names.difference(data.obs_names),]
	idx = tmp.X.sum(axis = 1).A1 < value
	pvec = tmp.X[idx, ].sum(axis = 0).A1
	pvec /= pvec.sum()
	data.uns['background_probs'] = pvec
	
	data.obsm['raw_probs'] = np.zeros((data.shape[0], adt.shape[1] + 1))
	data.obsm['raw_probs'][:, adt.shape[1]] = 1.0

	pool = multiprocessing.Pool(n_threads)
	data.obsm['raw_probs'][idx_df, :] = pool.map(estimate_probs_wrapper, adt_small)
	# data.obsm['raw_probs'][idx_df, :] = np.apply_along_axis(estimate_probs, 1, adt_small, pvec)

	demux_type = np.full(data.shape[0], 'unknown', dtype = 'object')
	assignments = np.full(data.shape[0], '', dtype = 'object')

	idx = data.obsm['raw_probs'][:,pvec.size] < unknown
	tmp = data.obsm['raw_probs'][idx,]
	norm_probs = tmp[:,0:pvec.size] / (1.0 - tmp[:,pvec.size])[:,None]
	idx_arr = norm_probs >= 0.1

	values1 = []
	values2 = []
	for i in range(idx_arr.shape[0]):
		droplet_type, droplet_id = get_droplet_info(idx_arr[i,])
		values1.append(droplet_type)
		values2.append(droplet_id)

	demux_type[idx] = values1
	data.obs['demux_type'] = pd.Categorical(demux_type, categories = ['singlet', 'doublet', 'unknown'])
	assignments[idx] = values2
	data.obs['assignment'] = pd.Categorical(assignments, categories = natsorted(np.unique(assignments)))
	
	end = time.time()
	print(end - start)
