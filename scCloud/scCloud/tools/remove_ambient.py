import numpy as np
import pandas as pd
import time
from natsort import natsorted

from scipy.sparse import vstack, hstack, csr_matrix
from scipy.stats import pearsonr

import warnings
warnings.filterwarnings("error")

def estimate_probs(arr, pvec, alpha_0, prob_noise, n_iter = 50):
	noise = pvec.size
	probs = np.zeros(pvec.size + 1)
	probs[-1] = prob_noise
	probs[:-1] = (1.0 - prob_noise) / pvec.size
	print(probs)
	# Beta distribution prior
	alpha = alpha_0 * prob_noise
	# EM algorithm
	z = np.zeros(pvec.size + 1)
	for i in range(n_iter):
		# E step
		z[:] = 0.0
		for j in range(pvec.size):
			if arr[j] > 0:
				p = probs[j] / (probs[noise] * pvec[j] + probs[j])
				z[j] += arr[j] * p
				z[noise] += arr[j] * (1.0 - p)
		# M step
		denom = z[:-1].sum()
		probs[noise] = (z[noise] + alpha) / (z[noise] + denom + alpha_0) # estiamte probs[noise]
		probs[:-1] = (1.0 - probs[noise]) * z[:-1] / denom
		print(probs)

	return probs

def estimate_probs2(arr, pvec, alpha = 0.5, alpha_noise = 1.0, n_iter = 50):
	idx_ori = arr > 0.0
	size = sum(idx_ori)
	arr = arr[idx_ori]
	pvec = pvec[idx_ori] / pvec[idx_ori].sum()

	probs = np.zeros(size + 1)
	z = np.zeros(size + 1)
	noise = size
	# set up alpha_vec
	alpha_vec = np.zeros(size + 1)
	alpha_vec[-1] = alpha_noise - 1.0
	alpha_vec[:-1][pvec > 0.0] = alpha - 1.0

	probs_mle = arr / arr.sum()
	probs[noise] = (probs_mle / (pvec + 1e-4)).min() + 0.01
	probs[:-1] = np.maximum(probs_mle - probs[noise] * pvec, 0.001)
	probs = probs / probs.sum()

	print(probs)
	# EM algorithm
	for i in range(n_iter):
		# E step
		z = np.copy(alpha_vec)
		for j in range(pvec.size):
			if arr[j] > 0:
				p = probs[j] / (probs[noise] * pvec[j] + probs[j])
				z[j] += arr[j] * p
				z[noise] += arr[j] * (1.0 - p)
		# M step
		old_probs = np.copy(probs)

		idx = z > 0.0
		probs[idx] = z[idx] / z[idx].sum()
		probs[~idx] = 0.0


		print("Diff = {}".format(np.linalg.norm(probs - old_probs, ord = 1)))
		print(probs)
	
	results = np.zeros(idx_ori.size + 1)
	results[-1] = probs[-1]
	results[:-1][idx_ori] = probs[:-1]
	return results
