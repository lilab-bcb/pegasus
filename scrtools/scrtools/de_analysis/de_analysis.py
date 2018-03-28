#!/usr/bin/env python

import numpy as np
from scipy.sparse import issparse
import fisher
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdr
import warnings
import pandas as pd
import scipy.stats as ss
import xlsxwriter

def collect_contingency_table(adata):
	clusts = adata.obs['louvain_groups'].value_counts()
	ct = np.zeros((adata.var_names.size, clusts.size, 2), dtype = np.uint)
	for i in range(clusts.size):
		label = clusts.index[i]
		count = clusts[i]
		mask = np.isin(adata.obs['louvain_groups'], label)
		ct[:, i, 0] = adata.X[mask].astype(bool).sum(axis = 0).A1
		ct[:, i, 1] = count - ct[:, i, 0]
	adata.uns["contingency_table"] = ct

def fisher_test(adata, clusts = None):
	if "contingency_table" not in adata.uns:
		collect_contingency_table(adata)
		print("Contingency table is collected.")

	ct = adata.uns["contingency_table"]
	token = ""

	if clusts != None:
		ct = ct[:, clusts, :]
		token = "_{0}".format(str(clusts))
	else:
		clusts = list(range(adata.obs['louvain_groups'].cat.categories.size))

	nclust = len(clusts)
	total = ct.sum(axis = 1)

	dtypes = [('percentage', np.dtype('float64')), ('fold_change', np.dtype('float64')), ('pval', np.dtype('float64')), ('qval', np.dtype('float64'))]
	
	for i in range(nclust):
		cpt = total - ct[:, i, :]
		pvals = fisher.pvalue_npy(ct[:, i, 0], ct[:, i, 1], cpt[:, 0], cpt[:, 1])[2]
		percents = ct[:, i, 0] / (ct[:, i, 0] + ct[:, i, 1]) * 100.0

		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			fold_change = ct[:, i, 0] / (ct[:, i, 0] + ct[:, i, 1]) / (cpt[:, 0] / (cpt[:, 0] + cpt[:, 1]))
			fold_change[np.isnan(fold_change)] = 0.0

		passed, qvals = fdr(pvals, alpha = 0.05)

		grt_idx = np.logical_and(passed, fold_change > 1.0)
		adata.uns["de_fisher{0}_{1}_up_genes".format(token, clusts[i])] = adata.var_names[grt_idx]
		values = list(zip(percents[grt_idx], fold_change[grt_idx], pvals[grt_idx], qvals[grt_idx]))
		adata.uns["de_fisher{0}_{1}_up_stats".format(token, clusts[i])] = np.array(values, dtype = dtypes)

		lsr_idx = np.logical_and(passed, fold_change < 1.0)
		adata.uns["de_fisher{0}_{1}_down_genes".format(token, clusts[i])] = adata.var_names[lsr_idx]
		values = list(zip(percents[lsr_idx], fold_change[lsr_idx], pvals[lsr_idx], qvals[lsr_idx]))
		adata.uns["de_fisher{0}_{1}_down_stats".format(token, clusts[i])] = np.array(values, dtype = dtypes)

		print("Cluster {0} is processed.".format(clusts[i]))



def t_test(adata, clusts = None):
	token = ""
	dtypes = [('mean_log_expression', np.dtype('float64')), ('log_fold_change', np.dtype('float64')), ('pval', np.dtype('float64')), ('qval', np.dtype('float64'))]

	if clusts != None:
		clusts = [str(x) for x in clusts] # convert to string
		token = "_{0}".format(str(clusts))
	else:
		clusts = adata.obs['louvain_groups'].cat.categories.tolist()

	nclust = len(clusts)
	
	n = n1 = n2 = 0

	mask = np.isin(adata.obs['louvain_groups'], clusts)
	mat = adata.X[mask]
	n = mat.shape[0]
	v = n - 2 # degree of freedom
	sm1 = mat.sum(axis = 0).A1 # sum of moment 1
	mat.data **= 2
	sm2 = mat.sum(axis = 0).A1 # sum of moment 2

	tvar = ss.t(v)
	for clust_id in clusts:
		mask = np.isin(adata.obs['louvain_groups'], clust_id)
		mat = adata.X[mask]
		n1 = mat.shape[0]
		n2 = n - n1
		sm1_1 = mat.sum(axis = 0).A1
		
		mean1 = sm1_1 / n1
		mean2 = (sm1 - sm1_1) / n2
		sp = ((sm2 - n1 * (mean1 ** 2) - n2 * (mean2 ** 2)) / v) ** 0.5

		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			tscore = (mean1 - mean2) / sp / ((1 / n1 + 1 / n2) ** 0.5)
			tscore[np.isnan(tscore)] = 1e3
			log_fold_change = mean1 / mean2
			log_fold_change[np.isnan(log_fold_change)] = 0.0

		pvals = tvar.sf(np.fabs(tscore)) * 2.0
		passed, qvals = fdr(pvals, alpha = 0.05)

		grt_idx = np.logical_and(passed, log_fold_change > 1.0)
		adata.uns["de_t{0}_{1}_up_genes".format(token, clust_id)] = adata.var_names[grt_idx]
		values = list(zip(mean1[grt_idx], log_fold_change[grt_idx], pvals[grt_idx], qvals[grt_idx]))
		adata.uns["de_t{0}_{1}_up_stats".format(token, clust_id)] = np.array(values, dtype = dtypes)

		lsr_idx = np.logical_and(passed, log_fold_change < 1.0)
		adata.uns["de_t{0}_{1}_down_genes".format(token, clust_id)] = adata.var_names[lsr_idx]
		values = list(zip(mean1[lsr_idx], log_fold_change[lsr_idx], pvals[lsr_idx], qvals[lsr_idx]))
		adata.uns["de_t{0}_{1}_down_stats".format(token, clust_id)] = np.array(values, dtype = dtypes)

		print("Cluster {0} is processed.".format(clust_id))


def write_results_to_excel(output_file, adata, test, clusts = None, threshold = 2.0):
	token = ""
	if clusts != None:
		token = "_{0}".format(str(clusts))
	else:
		clusts = list(range(adata.obs['louvain_groups'].cat.categories.size))

	thre_kw = 'fold_change' if test == 'fisher' else 'log_fold_change'
	sort_kw = 'percentage' if test == 'fisher' else 'mean_log_expression'

	writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
	for clust_id in clusts:
		df = pd.DataFrame(adata.uns["de_{0}{1}_{2}_up_stats".format(test, token, clust_id)], 
			index = pd.Index(adata.uns["de_{0}{1}_{2}_up_genes".format(test, token, clust_id)], name = "gene"))
		df = df.loc[df[thre_kw] >= threshold]
		df.sort_values(by = sort_kw, ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Cluster {0}, up-regulated".format(clust_id + 1))

		df = pd.DataFrame(adata.uns["de_{0}{1}_{2}_down_stats".format(test, token, clust_id)],
			index = pd.Index(adata.uns["de_{0}{1}_{2}_down_genes".format(test, token, clust_id)], name = "gene"))
		df = df.loc[df[thre_kw] <= 1.0 / threshold]
		df.sort_values(by = sort_kw, ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Cluster {0}, down-regulated".format(clust_id + 1))
	writer.save()
