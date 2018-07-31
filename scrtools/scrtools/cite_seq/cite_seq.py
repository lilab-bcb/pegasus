import numpy as np
import pandas as pd
import time
from natsort import natsorted

def estimate_probs(arr, pvec, n_iter = 200):
	probs = np.zeros(pvec.size + 1)
	z = np.zeros(pvec.size + 1)
	noise = pvec.size
	probs[:] = 0.1 / pvec.size
	probs[-1] = 0.9
	for i in range(n_iter):
		# E step
		z[:] = 0.0
		for j in range(pvec.size):
			if arr[j] > 0:
				p = probs[j] / (probs[noise] * pvec[j] + probs[j])
				z[j] += arr[j] * p
				z[noise] += arr[j] * (1.0 - p)
		# M step
		probs = z / z.sum()
	return probs

def get_droplet_type(arr):
	res = arr.nonzero()[0]
	return str(res[0] + 1) if res.size == 1 else 'doublet'

def demultiplex(data, adt):
	start = time.time()
	adt_small = adt[data.obs_names,].X.toarray()
	value = np.median(adt_small.sum(axis = 1))
	tmp = adt[adt.obs_names.difference(data.obs_names),]
	idx = tmp.X.sum(axis = 1).A1 < value
	pvec = tmp.X[idx, ].sum(axis = 0).A1
	pvec /= pvec.sum()
	data.obsm['raw_probs'] = np.apply_along_axis(estimate_probs, 1, adt_small, pvec)
	assignments = np.full(data.shape[0], 'unknown', dtype = 'object')
	idx = data.obsm['raw_probs'][:,pvec.size] < 0.5
	tmp = data.obsm['raw_probs'][idx,]
	norm_probs = tmp[:,0:pvec.size] / (1.0 - tmp[:,pvec.size])[:,None]
	idx_arr = norm_probs >= 0.1
	values = []
	for i in range(idx_arr.shape[0]):
		values.append(get_droplet_type(idx_arr[i,]))
	assignments[idx] = values
	data.obs['assignment'] = pd.Categorical(assignments, categories = natsorted(np.unique(assignments)))
	data.obs.loc[idx, 'assignment'] = values
	end = time.time()
	print(end - start)
