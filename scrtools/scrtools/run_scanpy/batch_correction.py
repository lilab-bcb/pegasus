#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import pandas as pd
import scanpy.api as sc
from scipy.sparse import issparse
from scipy.stats.mstats import gmean
from sklearn.preprocessing import StandardScaler
from collections import defaultdict

import scanpy_utils

def update_batch_id(data, func):
	scanpy_utils.add_obs_column(data, 'batch', func if func is not None else "lambda x: '-'.join(x.split('-')[:-2])")

def update_factor_id(data, func):
	scanpy_utils.add_obs_column(data, 'factor', func if func is not None else "lambda x: 'factor0'")

def standard_scaler(data, max_value = None):
	if issparse(data.X):
		data.X = data.X.toarray()

	scaler = StandardScaler(copy = False)	

	batch_ids = np.unique(data.obs['batch'])
	for batch_id in batch_ids:
		idx = np.isin(data.obs['batch'], batch_id)
		mat = data.X[idx]
		scaler.fit(mat)
		scaler.transform(mat)
		data.X[idx] = mat

	if max_value is not None:
		data.X[data.X > max_value] = max_value

	print("batch_correction.standard_scale done.")

def regress_out(data):
	if issparse(data.X):
		data.X = data.X.toarray()

	m = data.X.shape[1]
	fmeans = defaultdict(list)
	fstds = defaultdict(list)

	batch_ids = np.unique(data.obs['batch'])
	for batch_id in batch_ids:
		idx = np.isin(data.obs['batch'], batch_id)
		mat = data.X[idx]
		means = np.mean(mat, axis = 0)
		stds = np.std(mat, axis = 0)
		factor = data.obs['factor'][idx.nonzero()[0][0]]
		fmeans[factor].append(means)
		fstds[factor].append(stds)
		stds[stds < 1e-8] = 1.0 # avoid divide 0 error
		data.X[idx] = (mat - np.reshape(means, newshape = (1, m))) / np.reshape(stds, newshape = (1, m))

	factors = np.unique(data.obs['factor'])	
	for factor in factors:
		means = np.mean(fmeans[factor], axis = 0)
		stds = ma.getdata(gmean(ma.masked_less(fstds[factor], 1e-8), axis = 0))
		idx = np.isin(data.obs['factor'], factor)
		data.X[idx] =  data.X[idx] * np.reshape(stds, newshape = (1, m)) + np.reshape(means, newshape = (1, m))

	print("batch_correction.regress_out done.")

def filter_genes_dispersion(data, transformed, consider_batch, min_disp=0.5, max_disp=None, min_mean=0.0125, max_mean=7):	
	X = data.X.expm1() if transformed else data.X

	mean = X.mean(axis = 0).A1
	var = np.zeros(X.shape[1])

	if consider_batch:
		batch_ids = np.unique(data.obs['batch'])
		for batch_id in batch_ids:
			mat = X[np.isin(data.obs['batch'], batch_id)]
			m1 = mat.mean(axis = 0).A1
			m2 = mat.power(2).sum(axis = 0).A1
			var += m2 - mat.shape[0] * (m1 ** 2)
		factors = np.unique(data.obs['factor'])
		if factors.size > 1:
			for factor in factors:
				mat = X[np.isin(data.obs['factor'], factor)]
				var += mat.shape[0] * (mat.mean(axis = 0).A1 - mean) ** 2
		var /= X.shape[0] - 1 
	else:
		m2 = X.power(2).sum(axis = 0).A1
		var = (m2 - X.shape[0] * (mean ** 2)) / (X.shape[0] - 1)

	# now actually compute the dispersion
	mean[mean == 0] = 1e-12  # set entries equal to zero to small value
	dispersion = var / mean

	dispersion[dispersion == 0] = np.nan
	dispersion = np.log(dispersion)
	mean = np.log1p(mean)
	# all of the following quantities are "per-gene" here

	df = pd.DataFrame()
	df['mean'] = mean
	df['dispersion'] = dispersion
	df['mean_bin'] = pd.cut(df['mean'], bins=20)
	disp_grouped = df.groupby('mean_bin')['dispersion']
	disp_mean_bin = disp_grouped.mean()
	disp_std_bin = disp_grouped.std(ddof=1)
	df['dispersion_norm'] = (df['dispersion'].values  # use values here as index differs
							 - disp_mean_bin[df['mean_bin']].values) \
							 / disp_std_bin[df['mean_bin']].values
	dispersion_norm = df['dispersion_norm'].values.astype('float32')
	max_disp = np.inf if max_disp is None else max_disp
	dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
	gene_subset = np.logical_and.reduce((mean > min_mean, mean < max_mean,
										 dispersion_norm > min_disp,
										 dispersion_norm < max_disp))

	print("batch_correction.filter_genes_dispersion done.")
	print("{0} genes are selected.".format(gene_subset.sum()))

	return np.rec.fromarrays((gene_subset,
							  df['mean'].values,
							  df['dispersion'].values,
							  df['dispersion_norm'].values.astype('float32', copy=False)),
							  dtype=[('gene_subset', bool),
									 ('means', 'float32'),
									 ('dispersions', 'float32'),
									 ('dispersions_norm', 'float32')])
