import numpy as np
import pandas as pd
import time
from natsort import natsorted

import anndata
from scipy.sparse import vstack, hstack, csr_matrix
from scipy.stats import pearsonr

import warnings
warnings.filterwarnings("error")

import xlsxwriter

def estimate_probs_old(arr, pvec, n_iter = 200):
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

def estimate_probs(arr, pvec, alpha = 0.0, alpha_noise = 1.0, n_iter = 50):
	probs = np.zeros(pvec.size + 1)
	z = np.zeros(pvec.size + 1)
	noise = pvec.size
	# Estimate MLE without Generalized Dirichlet prior
	probs_mle = arr / arr.sum()
	probs[noise] = (probs_mle / pvec).min() + 0.01
	probs[:-1] = np.maximum(probs_mle - probs[noise] * pvec, 0.01)
	probs = probs / probs.sum()
	# EM algorithm
	for i in range(n_iter):
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
	# return results
	return probs

def get_droplet_type(arr):
	res = arr.nonzero()[0]
	return str(res[0] + 1) if res.size == 1 else 'doublet'

def demultiplex(data, adt, unknown_threshold = 0.5):
	start = time.time()

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
	
	data.obsm['raw_probs'] = np.zeros((data.shape[0], adt.shape[1] + 1))
	data.obsm['raw_probs'][:, adt.shape[1]] = 1.0
	data.obsm['raw_probs'][idx_df, :] = np.apply_along_axis(estimate_probs, 1, adt_small, pvec)

	assignments = np.full(data.shape[0], 'unknown', dtype = 'object')
	idx = data.obsm['raw_probs'][:,pvec.size] < unknown_threshold
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



def concatenate_ADTs(ADTs, names):
	obs_names = []
	Xs = []
	for adt, name in zip(ADTs, names):
		obs_names.append([name + '-' + x for x in adt.obs_names])
		Xs.append(adt.X)
	data = anndata.AnnData(X = vstack(Xs), obs = {"obs_names" : np.concatenate(obs_names)}, var = {"var_names" : ADTs[0].var_names})
	return data

def append_ADT_to_RNA_data(adt_data, df_a2c, data):
	adt_expr = np.zeros((data.shape[0], df_a2c.shape[0]))
	new_index, indexer = adt_data.obs_names.reindex(data.obs_names)
	idx = indexer >= 0
	adt_subset = adt_data[indexer[idx],]
	for i, row in df_a2c.iterrows():
		adt_expr[:, i] = np.maximum(np.log(adt_subset[:, row['antibody']].X + 1.0) - np.log(adt_subset[:, row['control']].X + 1.0), 0.0)
	
	new_data = anndata.AnnData(X = hstack([data.X, csr_matrix(adt_expr)], format = 'csr'), 
		obs = data.obs,
		obsm = data.obsm, 
		var = {"var_names" : np.concatenate([data.var_names, df_a2c['antibody'].apply(lambda x: 'AB-' + x)])})

	return new_data

def evaluate_concentration(data, df_r2a, attr, out_file):
	attrs = data.obs[attr].unique()
	results = []

	for rid, row in df_r2a.iterrows():
		genes = row['gene'].split(',')
		if len(genes) == 1:
			expr = data[:,genes].X.toarray()
		else:
			expr = np.log1p(np.expm1(data[:,genes].X.toarray()).sum(axis = 1))
		adts = data[:,row['antibody']].X.toarray()

		my_dict = {}
		for value in attrs:
			idx = np.isin(data.obs[attr], value)
			try:
				my_dict[value + '_rho'] = "{:.2f}".format(pearsonr(expr[idx], adts[idx])[0])
			except RuntimeWarning:
				my_dict[value + '_rho'] = "0.0"
			my_dict[value + '_mag_95'] = "{:.2f}".format(np.percentile(adts[idx], 95))
		results.append(my_dict)

	df = pd.DataFrame(results, index = df_r2a['antibody'])

	writer = pd.ExcelWriter(out_file, engine='xlsxwriter')
	df.to_excel(writer, sheet_name = "Eval")
	writer.save()
