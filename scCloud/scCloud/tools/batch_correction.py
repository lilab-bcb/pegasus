import time
import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy.sparse import issparse, csr_matrix
from collections import defaultdict

from numba import njit
from math import sqrt

import skmisc.loess as sl
from scCloud.plotting import plot_hvg

def set_group_attribute(data, attribute_string):
	if attribute_string is None:
		data.obs['Group'] = 'one_group'
	elif attribute_string.find('=') >= 0:
		attr, value_str = attribute_string.split('=')
		assert attr in data.obs.columns
		values = value_str.split(';')
		data.obs['Group'] = '0'
		for group_id, value in enumerate(values):
			vals = value.split(',')
			idx = np.isin(data.obs[attr], vals)
			data.obs.loc[idx, 'Group'] = str(group_id + 1)
	elif attribute_string.find('+') >= 0:
		attrs = attribute_string.split('+')
		assert np.isin(attrs, data.obs.columns).sum() == len(attrs)
		data.obs['Group'] = data.obs[attrs].apply(lambda x: '+'.join(x), axis = 1)
	else:
		assert attribute_string in data.obs.columns
		data.obs['Group'] = data.obs[attribute_string]



def collect_highly_variable_gene_matrix(data):
	data_c = data[:, data.var['highly_variable_genes']].copy()
	data_c.X = data_c.X.toarray()
	return data_c



def estimate_adjustment_matrices(data):
	start = time.time()

	channels = data.obs['Channel'].unique()
	if channels.size <= 1:
		print("Warning: data only contains 1 channel. Batch correction disabled!")
		return None

	means = np.zeros((data.shape[1], channels.size))
	partial_sum = np.zeros((data.shape[1], channels.size))
	ncells = np.zeros(channels.size)
	stds = np.zeros((data.shape[1], channels.size))

	group_dict = defaultdict(list)
	assert issparse(data.X)

	for i, channel in enumerate(channels):
		idx = np.isin(data.obs['Channel'], channel)
		mat = data.X[idx].astype(np.float64)

		ncells[i] = mat.shape[0]
		if ncells[i] == 1:
			means[:, i] = mat.toarray()[0]
		else:
			means[:, i] = mat.mean(axis = 0).A1
			m2 = mat.power(2).sum(axis = 0).A1
			partial_sum[:, i] = m2 - ncells[i] * (means[:, i] ** 2)

		group = data.obs['Group'][idx.nonzero()[0][0]]
		group_dict[group].append(i)

	partial_sum[partial_sum < 1e-6] = 0.0

	overall_means = np.dot(means, ncells) / data.shape[0]
	batch_adjusted_vars = np.zeros(data.shape[1])
	groups = data.obs['Group'].unique()
	for group in groups:
		gchannels = group_dict[group]
		ncell = ncells[gchannels].sum()
		gm = np.dot(means[:, gchannels], ncells[gchannels]) / ncell
		gs = partial_sum[:, gchannels].sum(axis = 1) / ncell

		if groups.size > 1:
			batch_adjusted_vars += ncell * ((gm - overall_means) ** 2)

		for i in gchannels:
			if ncells[i] > 1:
				stds[:, i] = (partial_sum[:, i] / (ncells[i] - 1.0)) ** 0.5
			outliers = stds[:, i] < 1e-6
			normals = np.logical_not(outliers)
			stds[outliers, i] = 1.0
			stds[normals, i] = gs[normals] / stds[normals, i]
			means[:, i] = gm - stds[:, i] * means[:, i]

	data.uns['Channels'] = channels
	data.uns['Groups'] = groups
	data.varm['means'] = means
	data.varm['stds'] = stds

	data.var['ba_mean'] = overall_means
	data.var['ba_var'] = (batch_adjusted_vars + partial_sum.sum(axis = 1)) / (data.shape[0] - 1.0)

	end = time.time()
	print("batch_correction.estimate_adjustment_matrices is done. Time spent = {:.2f}s.".format(end - start))



def correct_batch_effects(data):
	start = time.time()

	if 'Channels' not in data.uns:
		print("Warning: Batch correction parameters are not calculated. Batch correction skipped!")
		return None

	assert not issparse(data.X)
	m = data.shape[1]
	for i, channel in enumerate(data.uns['Channels']):
		idx = np.isin(data.obs['Channel'], channel)
		if idx.sum() == 0:
			continue
		data.X[idx] = data.X[idx] * np.reshape(data.varm['stds'][:, i], newshape = (1, m)) + np.reshape(data.varm['means'][:, i], newshape = (1, m))
	data.X[data.X < 0.0] = 0.0

	end = time.time()
	print("batch_correction.correct_batch_effects done. Time spent = {:.2f}s".format(end - start))



def select_highly_variable_genes(data, consider_batch, flavor = 'scCloud', n_top = 2000, span = 0.02, min_disp = 0.5, max_disp = np.inf, min_mean = 0.0125, max_mean = 7, plot_hvg_fig = None):
	start = time.time()

	if consider_batch and 'Channels' not in data.uns:
		print("Warning: Batch correction parameters are not calculated. Switch to not considering batch for variable gene selection.")
		consider_batch = False


	robust_idx = data.var['robust'].values
	hvg_index = np.zeros(robust_idx.sum(), dtype = bool)

	if flavor == 'scCloud':
		if consider_batch:
			mean = data.var.loc[robust_idx, 'ba_mean']
			var = data.var.loc[robust_idx, 'ba_var']
		else:
			X = data.X[:, robust_idx]
			mean = X.mean(axis = 0).A1
			m2 = X.power(2).sum(axis = 0).A1	
			var = (m2 - X.shape[0] * (mean ** 2)) / (X.shape[0] - 1)

		lobj = sl.loess(mean, var, span = span, degree = 2)
		lobj.fit()
		
		rank1 = np.zeros(hvg_index.size, dtype = int)
		rank2 = np.zeros(hvg_index.size, dtype = int)

		delta = var - lobj.outputs.fitted_values
		fc = var / lobj.outputs.fitted_values

		rank1[np.argsort(delta)[::-1]] = range(hvg_index.size)
		rank2[np.argsort(fc)[::-1]] = range(hvg_index.size)
		rank = rank1 + rank2

		data.var['hvg_rank'] = -1
		data.var.loc[robust_idx, 'hvg_rank'] = rank

		hvg_index[np.argsort(rank)[:n_top]] = True

		if plot_hvg_fig is not None:
			plot_hvg(mean, var, lobj.outputs.fitted_values, hvg_index, plot_hvg_fig + '.hvg.pdf')

	else:
		assert flavor == 'Seurat'
		X = data.X[:, robust_idx].expm1()
		mean = X.mean(axis = 0).A1
		var = np.zeros(X.shape[1])

		if consider_batch:
			for channel in data.uns['Channels']:
				mat = X[np.isin(data.obs['Channel'], channel)]
				if mat.shape[0] == 0:
					continue
				m1 = mat.mean(axis = 0).A1
				m2 = mat.power(2).sum(axis = 0).A1
				var += m2 - mat.shape[0] * (m1 ** 2)

			if data.uns['Groups'].size > 1:
				for group in data.uns['Groups']:
					mat = X[np.isin(data.obs['Group'], group)]
					if mat.shape[0] == 0:
						continue
					var += mat.shape[0] * ((mat.mean(axis = 0).A1 - mean) ** 2)
			var /= X.shape[0] - 1.0 
		else:
			m2 = X.power(2).sum(axis = 0).A1
			var = (m2 - X.shape[0] * (mean ** 2)) / (X.shape[0] - 1)

		dispersion = np.full(X.shape[1], np.nan)
		idx_valid = (mean > 0.0) & (var > 0.0)
		dispersion[idx_valid] = var[idx_valid] / mean[idx_valid]

		mean = np.log1p(mean)
		dispersion = np.log(dispersion)

		df = pd.DataFrame({'log_dispersion' : dispersion, 'bin' : pd.cut(mean, bins = 20)})
		log_disp_groups = df.groupby('bin')['log_dispersion']
		log_disp_mean = log_disp_groups.mean()
		log_disp_std = log_disp_groups.std(ddof = 1)
		log_disp_zscore = (df['log_dispersion'].values - log_disp_mean.loc[df['bin']].values) / log_disp_std.loc[df['bin']].values
		log_disp_zscore[np.isnan(log_disp_zscore)] = 0.0

		hvg_index = np.logical_and.reduce((mean > min_mean, mean < max_mean, log_disp_zscore > min_disp, log_disp_zscore < max_disp))
	
	data.var.loc[robust_idx, 'highly_variable_genes'] = hvg_index

	end = time.time()
	print("batch_correction.select_highly_variable_genes done. Time spent = {:.2f}s.".format(end - start))
	print("{0} genes are selected.".format(data.var['highly_variable_genes'].sum()))
