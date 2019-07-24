#!/usr/bin/env python

import time
import numpy as np
import pandas as pd
from natsort import natsorted
from scipy.sparse import csc_matrix
import scipy.stats as ss
from scipy.stats import f_oneway
from statsmodels.stats.multitest import fdrcorrection as fdr
from joblib import Parallel, delayed
from scCloud.tools import read_input



def search_genes(data, gene_list, measure = 'percentage'):
	"""Extract and display gene expressions for each cluster from an `anndata` object.

	This function helps to see marker expressions in clusters via the interactive python environment.

	Parameters
	----------

	data : `anndata` object
		An `anndata` object containing the expression matrix and differential expression results.
	gene_list : `list[str]`
		A list of gene symbols.
	measure : `str`
		Can be either `percentage` or `mean_log_expression`. `percentage` shows the percentage of cells expressed the genes and `mean_log_expression` shows the mean log expression.

	Returns
	-------
	`pandas.DataFrame`
		A data frame containing marker expressions in each cluster.

	Examples
	--------
	>>> results = misc.search_genes(data, ['CD3E', 'CD4', 'CD8'], measure = 'percentage')
	"""

	ncluster =  sum([1 for x in data.var.columns if x.startswith('mean_log_expression_')])
	columns = ["{}_{}".format(measure, x + 1) for x in range(ncluster)]
	return data.var.reindex(index = gene_list, columns = columns)

def search_de_genes(data, gene_list, test = 'fisher', thre = 1.5):
	"""Extract and display differential expression analysis results of markers for each cluster from an `anndata` object.

	This function helps to see if markers are up or down regulated in each cluster via the interactive python environment. `++` indicates up-regulated and fold change >= threshold, `+` indicates up-regulated but fold change < threshold, `--` indicates down-regulated and fold change <= 1 / threshold, `-` indicates down-regulated but fold change > 1 / threshold, '?' indicates not differentially expressed.

	Parameters
	----------

	data : `anndata` object
		An `anndata` object containing the expression matrix and differential expression results.
	gene_list : `list[str]`
		A list of gene symbols.
	test : `str`, optional (default: `fisher`)
		Differential expression test to look at, could be either `t`, `fisher` or `mwu`.
	thre : `float`, optional (default: `1.5`)
		Fold change threshold to determine if the marker is a strong DE (`++` or `--`) or weak DE (`+` or `-`).

	Returns
	-------
	`pandas.DataFrame`
		A data frame containing marker differential expression results for each cluster.

	Examples
	--------
	>>> results = misc.search_de_genes(data, ['CD3E', 'CD4', 'CD8'], test = 'fisher', thre = 2.0)
	"""

	ngene = len(gene_list)
	ncluster = sum([1 for x in data.var.columns if x.startswith('mean_log_expression_')])
	results = np.zeros((ngene, ncluster), dtype = np.dtype('U4'))
	columns = [str(x + 1) for x in range(ncluster)]
	df_de = data.var.reindex(index = gene_list, columns = [test + "_qval_" + x for x in columns])
	if test == 'fisher':
		df_fc = data.var.reindex(index = gene_list, columns = ["percentage_fold_change_" + x for x in columns])
	else:
		df_fc = np.exp(data.var.reindex(index = gene_list, columns = ["log_fold_change_" + x for x in columns]))
	results[:] = '?'
	results[np.isnan(df_de)] = 'NaN'
	results[(df_de <= 0.05).values & (df_fc > 1.0).values] = '+'
	results[(df_de <= 0.05).values & (df_fc >= thre).values] = '++'
	results[(df_de <= 0.05).values & (df_fc < 1.0).values] = '-'
	results[(df_de <= 0.05).values & (df_fc <= 1.0 / thre).values] = '--'
	df = pd.DataFrame(results, index = gene_list, columns = columns)
	return df

def show_attributes(input_file, show_attributes, show_gene_attributes, show_values_for_attributes):
	data = read_input(input_file, mode = 'r')
	if show_attributes:
		print("Available sample attributes in input dataset: {0}".format(', '.join(data.obs.columns.values)))
	if show_gene_attributes:
		print("Available gene attributes in input dataset: {0}".format(', '.join(data.var.columns.values)))
	if not show_values_for_attributes is None:
		for attr in show_values_for_attributes.split(','):
			print("Available values for attribute {0}: {1}.".format(attr, ', '.join(np.unique(data.obs[attr]))))



def perform_oneway_anova(data, glist, restriction_vec, group_str, fdr_alpha = 0.05):
	selected = np.ones(data.shape[0], dtype = bool)
	for rest_str in restriction_vec:
		attr, value_str = rest_str.split(':')
		values = value_str.split(',')
		selected = selected & np.isin(data.obs[attr], values)
	gene_list = np.array(glist)
	gene_list = gene_list[np.isin(gene_list, data.var_names)]
	newdat = data[selected, :][:, gene_list].copy()
	newdat.X = newdat.X.toarray()
	group_attr, tmp_str = group_str.split(':')
	groups_str = tmp_str.split(';')
	ngr = len(groups_str)
	group_names = []
	group_idx = np.zeros((ngr, newdat.shape[0]), dtype = bool)
	for i, gstr in enumerate(groups_str):
		name, values = gstr.split('~')
		group_names.extend([name + '_mean', name + '_percent'])
		group_idx[i] = np.isin(newdat.obs[group_attr], values.split(','))
	np.warnings.filterwarnings('ignore')
	stats = np.zeros((len(gene_list), 3 + ngr * 2))
	for i in range(len(gene_list)):
		arr_list = []
		for j in range(group_idx.shape[0]):
			arr = newdat.X[group_idx[j], i]
			stats[i, 3 + j * 2] = arr.mean()
			stats[i, 3 + j * 2 + 1] = (arr > 0).sum() * 100.0 / arr.size
			arr_list.append(arr)
		stats[i, 0], stats[i, 1] = f_oneway(*arr_list)
		if np.isnan(stats[i, 0]):
			stats[i, 0] = 0.0
			stats[i, 1] = 1.0
	passed, stats[:, 2] = fdr(stats[:, 1])
	cols = ['fstat', 'pval', 'qval']
	cols.extend(group_names)
	raw_results = pd.DataFrame(stats, columns = cols, index = gene_list)
	results = raw_results[raw_results['qval'] <= fdr_alpha]
	results = results.sort_values('qval')
	return results, raw_results



# labels, cluster labels for each sample; conds, conditions; cond_order, condition orders
def calc_mwu(clust_label, labels, conds, cond_order, gene_names, data, indices, indptr, shape):
	csc_mat = csc_matrix((data, indices, indptr), shape = shape)
	ngene = shape[1]
	log_fc = np.zeros(ngene)
	U_stats = np.zeros(ngene)
	pvals = np.zeros(ngene)

	idx = labels == clust_label
	exprs = np.zeros(idx.sum())

	idx_x = conds[idx] == cond_order[0]
	idx_y = conds[idx] == cond_order[1]

	local_mat = csc_mat[idx, :]

	for j in range(ngene):
		vec = local_mat[:, j]
		if vec.size > 0:
			exprs[vec.indices] = vec.data
			log_fc[j] = np.mean(exprs[idx_x]) - np.mean(exprs[idx_y])
			U_stats[j], pvals[j] = ss.mannwhitneyu(exprs[idx_x], exprs[idx_y], alternative = 'two-sided')
		else:
			log_fc[j] = 0.0
			U_stats[j] = 0.0
			pvals[j] = 1.0
		exprs[:] = 0.0

	passed, qvals = fdr(pvals)

	df = pd.DataFrame({"log_fc": log_fc,
					   "mwu_U": U_stats,
					   "mwu_pval": pvals,
					   "mwu_qval": qvals},
					   index = gene_names)

	print("Cluster {0} is processed.".format(clust_label))

	return df



def mwu_test(data, attr_label, attr_cond, cond_order = None, n_jobs = 1, temp_folder = None):
	start = time.time()

	csc_mat = data.X.tocsc()

	if cond_order is None:
		cond_order = data.obs[attr_cond].cat.categories.values

	results = Parallel(n_jobs = n_jobs, max_nbytes = 1e7, temp_folder = temp_folder)(delayed(calc_mwu)(clust_label, data.obs[attr_label].values, data.obs[attr_cond].values, cond_order, data.var_names.values, csc_mat.data, csc_mat.indices, csc_mat.indptr, csc_mat.shape) for clust_label in data.obs[attr_label].cat.categories)
	result_dict = {x : y for x, y in zip(data.obs[attr_label].cat.categories, results)}

	end = time.time()
	print("Mann-Whitney U test is done. Time spent = {:.2f}s.".format(end - start))

	return result_dict
