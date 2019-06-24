import time
import numpy as np
from scipy.sparse import issparse
from collections import defaultdict



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
