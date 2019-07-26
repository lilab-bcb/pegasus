import time
import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix, csc_matrix
import scipy.stats as ss
import fisher
from statsmodels.stats.multitest import fdrcorrection as fdr
import sklearn.metrics as sm
import warnings
import xlsxwriter
from natsort import natsorted
from joblib import Parallel, delayed

from . import read_input, non_de_attrs


def collect_contingency_table(data, X, labels='louvain_labels'):
	if data.obs[labels].dtype.name != 'category':
		print('Warning: {} is not categorical, scCloud will convert it to categorical.'.format(labels))
		data.obs[labels] = pd.Categorical(data.obs[labels])

	nclust = data.obs[labels].cat.categories.size
	data.uns["de_nclust"] = nclust
	ct = np.zeros((data.var_names.size, nclust, 2), dtype=np.uint)
	for i, label in enumerate(data.obs[labels].cat.categories):
		mask = np.isin(data.obs[labels], label)
		ct[:, i, 0] = X[mask].getnnz(axis=0)
		ct[:, i, 1] = mask.sum() - ct[:, i, 0]
	data.uns["contingency_table"] = ct


def calc_stat_and_t(i, clust_label, labels, gene_names, data, indices, indptr, shape, sm1, sm2, ct, total):
	mat = csr_matrix((data, indices, indptr), shape=shape)

	n = shape[0]
	pvals = np.zeros(shape[1])
	qvals = np.zeros(shape[1])
	percent_fold_change = np.zeros(shape[1])

	mask = labels == clust_label
	clust_mat = mat[mask]
	n1 = clust_mat.shape[0]
	n2 = n - n1

	assert n1 > 0 and n2 > 0

	pvals[:] = 1.01
	qvals[:] = 1.01

	sm1_1 = clust_mat.sum(axis=0).A1
	mean1 = sm1_1 / n1
	mean2 = (sm1 - sm1_1) / n2

	if n1 > 1 and n2 > 1:
		sm2_1 = clust_mat.power(2).sum(axis=0).A1
		s1sqr = (sm2_1 - n1 * (mean1 ** 2)) / (n1 - 1)
		s2sqr = ((sm2 - sm2_1) - n2 * (mean2 ** 2)) / (n2 - 1)

		var_est = s1sqr / n1 + s2sqr / n2

		idx = var_est > 0.0
		if idx.sum() > 0:
			tscore = (mean1[idx] - mean2[idx]) / np.sqrt(var_est[idx])
			v = (var_est[idx] ** 2) / ((s1sqr[idx] / n1) ** 2 / (n1 - 1) + (s2sqr[idx] / n2) ** 2 / (n2 - 1))
			pvals[idx] = ss.t.sf(np.fabs(tscore), v) * 2.0  # two-sided
			passed, qvals[idx] = fdr(pvals[idx])

	# calculate WAD, Weighted Average Difference, https://almob.biomedcentral.com/articles/10.1186/1748-7188-3-8
	log_fold_change = mean1 - mean2
	x_avg = (mean1 + mean2) / 2
	x_max = x_avg.max()
	x_min = x_avg.min() - 0.001  # to avoid divide by zero
	weights = (x_avg - x_min) / (x_max - x_min)
	wads = log_fold_change * weights

	# calculate percentage expressed and percent fold change
	percents = ct[:, i, 0] / (ct[:, i, 0] + ct[:, i, 1]) * 100.0
	cpt = total - ct[:, i, :]
	percents_other = cpt[:, 0] / (cpt[:, 0] + cpt[:, 1]) * 100.0

	idx = percents > 0.0
	idx_other = percents_other > 0.0
	percent_fold_change[(~idx) & (~idx_other)] = 0.0
	percent_fold_change[idx & (~idx_other)] = np.inf
	percent_fold_change[idx_other] = percents[idx_other] / percents_other[idx_other]

	df = pd.DataFrame({"percentage_{0}".format(clust_label): percents,
					   "percentage_other_{0}".format(clust_label): percents_other,
					   "mean_log_expression_{0}".format(clust_label): mean1,
					   "percentage_fold_change_{0}".format(clust_label): percent_fold_change,
					   "log_fold_change_{0}".format(clust_label): log_fold_change,
					   "WAD_score_{0}".format(clust_label): wads,
					   "t_pval_{0}".format(clust_label): pvals,
					   "t_qval_{0}".format(clust_label): qvals},
					  index=gene_names)

	print("Cluster {0} is processed.".format(clust_label))

	return df


def collect_stat_and_t_test(data, X, clusts=None, labels='louvain_labels', n_jobs=1, temp_folder=None):
	start = time.time()

	collect_contingency_table(data, X, labels=labels)
	print("Contingency table is collected.")

	ct = data.uns["contingency_table"]

	if clusts is None:
		idx = (data.obs[labels].value_counts(sort=False) > 0).values
	else:
		idx = np.isin(data.obs[labels].cat.categories, clusts)

	if idx.sum() < data.uns["de_nclust"]:
		ct = ct[:, idx, :]

	clusts = data.obs[labels].cat.categories.values[idx]
	total = ct.sum(axis=1)

	mask = np.isin(data.obs[labels], clusts)
	mat = X[mask]
	sm1 = mat.sum(axis=0).A1  # sum of moment 1
	sm2 = mat.power(2).sum(axis=0).A1  # sum of moment 2

	result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
		delayed(calc_stat_and_t)(i, clusts[i], data.obs[labels][mask].values, data.var_names.values, mat.data,
								 mat.indices, mat.indptr, mat.shape, sm1, sm2, ct, total) for i in range(clusts.size))

	end = time.time()
	print("Welch's t-test is done. Time spent = {:.2f}s.".format(end - start))

	return result_list


def calc_fisher(i, clust_label, gene_names, ct, total):
	cpt = total - ct[:, i, :]
	pvals = fisher.pvalue_npy(ct[:, i, 0], ct[:, i, 1], cpt[:, 0], cpt[:, 1])[2]
	passed, qvals = fdr(pvals)
	df = pd.DataFrame({"fisher_pval_{0}".format(clust_label): pvals,
					   "fisher_qval_{0}".format(clust_label): qvals},
					  index=gene_names)

	print("Cluster {0} is processed.".format(clust_label))

	return df


# clusts is a list of strings
def fisher_test(data, X, clusts=None, labels='louvain_labels', n_jobs=1, temp_folder=None):
	start = time.time()

	if "contingency_table" not in data.uns:
		collect_contingency_table(data, X, labels=labels)
		print("Contingency table is collected.")

	ct = data.uns["contingency_table"]

	if clusts is None:
		idx = (data.obs[labels].value_counts(sort=False) > 0).values
	else:
		idx = np.isin(data.obs[labels].cat.categories, clusts)

	if idx.sum() < data.uns["de_nclust"]:
		ct = ct[:, idx, :]

	clusts = data.obs[labels].cat.categories.values[idx]
	total = ct.sum(axis=1)

	result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
		delayed(calc_fisher)(i, clusts[i], data.var_names.values, ct, total) for i in range(clusts.size))

	end = time.time()
	print("Fisher's exact test is done. Time spent = {:.2f}s.".format(end - start))

	return result_list


def calc_mwu(clust_label, labels, gene_names, data, indices, indptr, shape):
	csc_mat = csc_matrix((data, indices, indptr), shape=shape)
	nsample = shape[0]
	ngene = shape[1]

	idx_x = labels == clust_label
	idx_y = ~idx_x

	exprs = np.zeros(nsample)
	U_stats = np.zeros(ngene)
	pvals = np.zeros(ngene)

	for j in range(ngene):
		exprs[:] = 0.0
		vec = csc_mat[:, j]
		if vec.size > 0:
			exprs[vec.indices] = vec.data
			U_stats[j], pvals[j] = ss.mannwhitneyu(exprs[idx_x], exprs[idx_y], alternative='two-sided')
		else:
			U_stats[j] = 0.0
			pvals[j] = 1.0
	passed, qvals = fdr(pvals)

	df = pd.DataFrame({"mwu_U_{0}".format(clust_label): U_stats,
					   "mwu_pval_{0}".format(clust_label): pvals,
					   "mwu_qval_{0}".format(clust_label): qvals},
					  index=gene_names)

	print("Cluster {0} is processed.".format(clust_label))

	return df


def mwu_test(data, X, clusts=None, labels='louvain_labels', n_jobs=1, temp_folder=None):
	start = time.time()

	if clusts is None:
		idx = (data.obs[labels].value_counts(sort=False) > 0).values
	else:
		idx = np.isin(data.obs[labels].cat.categories, clusts)

	clusts = data.obs[labels].cat.categories.values[idx]

	mask = np.isin(data.obs[labels], clusts)
	csc_mat = X[mask].tocsc()

	result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
		delayed(calc_mwu)(clust_label, data.obs[labels][mask].values, data.var_names.values, csc_mat.data,
						  csc_mat.indices, csc_mat.indptr, csc_mat.shape) for clust_label in clusts)

	end = time.time()
	print("Mann-Whitney U test is done. Time spent = {:.2f}s.".format(end - start))

	return result_list


def calc_roc(clust_label, labels, gene_names, data, indices, indptr, shape):
	csc_mat = csc_matrix((data, indices, indptr), shape=shape)
	nsample = shape[0]
	ngene = shape[1]

	mask = labels == clust_label
	exprs = np.zeros(nsample)

	auc = np.zeros(ngene)
	predpower = np.zeros(ngene)
	tpr_at_fpr01 = np.zeros(ngene)
	tpr_at_fpr025 = np.zeros(ngene)
	tpr_at_fpr03 = np.zeros(ngene)
	tpr_at_fpr05 = np.zeros(ngene)

	for j in range(ngene):
		exprs[:] = 0.0
		vec = csc_mat[:, j]
		exprs[vec.indices] = vec.data
		fpr, tpr, thresholds = sm.roc_curve(mask, exprs)

		auc[j] = sm.auc(fpr, tpr)
		predpower[j] = np.fabs(auc[j] - 0.5) * 2

		tpr_at_fpr01[j] = tpr[np.argmax(fpr[fpr < 0.1])]
		tpr_at_fpr025[j] = tpr[np.argmax(fpr[fpr < 0.25])]
		tpr_at_fpr03[j] = tpr[np.argmax(fpr[fpr < 0.3])]
		tpr_at_fpr05[j] = tpr[np.argmax(fpr[fpr < 0.5])]

	df = pd.DataFrame({"auc_{0}".format(clust_label): auc,
					   "predpower_{0}".format(clust_label): predpower,
					   "tpr_at_fpr01_{0}".format(clust_label): tpr_at_fpr01,
					   "tpr_at_fpr025_{0}".format(clust_label): tpr_at_fpr025,
					   "tpr_at_fpr03_{0}".format(clust_label): tpr_at_fpr03,
					   "tpr_at_fpr05_{0}".format(clust_label): tpr_at_fpr05},
					  index=gene_names)

	print("Cluster {0} is processed.".format(clust_label))

	return df


def calc_roc_stats(data, X, clusts=None, labels='louvain_labels', n_jobs=1, temp_folder=None):
	start = time.time()

	if clusts is None:
		idx = (data.obs[labels].value_counts(sort=False) > 0).values
	else:
		idx = np.isin(data.obs[labels].cat.categories, clusts)

	clusts = data.obs[labels].cat.categories.values[idx]

	mask = np.isin(data.obs[labels], clusts)
	csc_mat = X[mask].tocsc()

	result_list = Parallel(n_jobs=n_jobs, max_nbytes=1e7, temp_folder=temp_folder)(
		delayed(calc_roc)(clust_label, data.obs[labels][mask].values, data.var_names.values, csc_mat.data,
						  csc_mat.indices, csc_mat.indptr, csc_mat.shape) for clust_label in clusts)

	end = time.time()
	print("ROC statistics are calculated. Time spent = {:.2f}s.".format(end - start))

	return result_list


def format_short_output_cols(df, cols_short_format):
	""" Round related float columns to 3 decimal points."""
	cols_short_format_idx = [df.columns.get_loc(c) for c in df.columns if c in cols_short_format]
	df.iloc[:, cols_short_format_idx] = df.iloc[:, cols_short_format_idx].round(3)
	return df


test2fields = {'t': ['t_pval', 't_qval'], 'fisher': ['fisher_pval', 'fisher_qval'],
			   'mwu': ['mwu_U', 'mwu_pval', 'mwu_qval']}


def write_results_to_excel(output_file, df, alpha=0.05):
	clusts = natsorted([x[10:] for x in df.columns if x.startswith("WAD_score_")])
	tests = [x for x in ['t', 'fisher', 'mwu'] if "{0}_qval_{1}".format(x, clusts[0]) in df.columns]
	has_roc = "auc_{0}".format(clusts[0]) in df.columns

	cols = ["percentage", "percentage_other", "percentage_fold_change", "mean_log_expression", "log_fold_change",
			"WAD_score"]
	if has_roc:
		cols.extend(["auc", "predpower"])
	if has_roc:
		cols.extend(["tpr_at_fpr01", "tpr_at_fpr025", "tpr_at_fpr03", "tpr_at_fpr05"])
	cols_short_format = cols.copy()
	for test in tests:
		cols.extend(test2fields[test])

	workbook = xlsxwriter.Workbook(output_file, {'nan_inf_to_errors': True})
	workbook.formats[0].set_font_size(9)
	for clust_id in clusts:
		idx = df["{0}_qval_{1}".format(tests[0], clust_id)] <= alpha
		for test in tests[1:]:
			idx = idx & (df["{0}_qval_{1}".format(test, clust_id)] <= alpha)

		idx_up = idx & (df["WAD_score_{0}".format(clust_id)] > 0.0)
		idx_down = idx & (df["WAD_score_{0}".format(clust_id)] < 0.0)
		assert idx_up.sum() + idx_down.sum() == idx.sum()

		col_names = ["{0}_{1}".format(x, clust_id) for x in cols]
		df_up = pd.DataFrame(df.loc[idx_up.values, col_names])
		df_up.rename(columns=lambda x: '_'.join(x.split('_')[:-1]), inplace=True)
		df_up.sort_values(by="auc" if has_roc else "WAD_score", ascending=False, inplace=True)
		# format output as excel table
		df_up = format_short_output_cols(df_up, cols_short_format)
		worksheet = workbook.add_worksheet(name="{0} up".format(clust_id))
		df_up.reset_index(inplace=True)
		df_up.rename(index=str, columns={"index": "gene"}, inplace=True)
		if len(df_up.index) > 0:
			worksheet.add_table(0, 0, len(df_up.index), len(df_up.columns) - 1,
								{'data': np.array(df_up), 'style': 'Table Style Light 1',
								 'first_column': True, 'header_row': True,
								 'columns': [{'header': x} for x in df_up.columns.values]})
		else:
			worksheet.write_row(0, 0, df_up.columns.values)

		df_down = pd.DataFrame(df.loc[idx_down.values, col_names])
		df_down.rename(columns=lambda x: '_'.join(x.split('_')[:-1]), inplace=True)
		df_down.sort_values(by="auc" if has_roc else "WAD_score", ascending=True, inplace=True)
		# format output as excel table
		worksheet = workbook.add_worksheet(name="{0} dn".format(clust_id))
		df_down = format_short_output_cols(df_down, cols_short_format)
		df_down.reset_index(inplace=True)
		df_down.rename(index=str, columns={"index": "gene"}, inplace=True)
		if len(df_up.index) > 0:
			worksheet.add_table(0, 0, len(df_down.index), len(df_down.columns) - 1,
								{'data': np.array(df_down), 'style': 'Table Style Light 1',
								 'first_column': True, 'header_row': True,
								 'columns': [{'header': x} for x in df_down.columns.values]})
		else:
			worksheet.write_row(0, 0, df_down.columns.values)
	workbook.close()

	print("Excel spreadsheet is written.")


def markers(data, labels, n_jobs, run_fisher, run_mwu, run_roc):
	__markers(data, data.X, labels, n_jobs=n_jobs, run_fisher=run_fisher, run_mwu=run_mwu, run_roc=run_roc)

def __markers(data, X, labels, n_jobs=1, run_fisher=True, run_mwu=True, run_roc=True, temp_folder=None):
	non_de = [x for x in non_de_attrs if x in data.var]
	de_results = [data.var[non_de]]

	print("Begin t_test.")
	de_results.extend(collect_stat_and_t_test(data, X, labels=labels, n_jobs=n_jobs, temp_folder=temp_folder))

	if run_fisher:
		print("Begin Fisher's exact test.")
		de_results.extend(fisher_test(data, X, labels=labels, n_jobs=n_jobs, temp_folder=temp_folder))

	if run_mwu:
		print("Begin Mann-Whitney U test.")
		de_results.extend(mwu_test(data, X, labels=labels, n_jobs=n_jobs, temp_folder=temp_folder))

	if run_roc:
		print("Begin calculating ROC statistics.")
		de_results.extend(calc_roc_stats(data, X, labels=labels, n_jobs=n_jobs, temp_folder=temp_folder))

	data.var = pd.concat(de_results, axis=1)
	data.uns['de_labels'] = labels

def run_de_analysis(input_file, output_excel_file, labels, n_jobs, alpha, run_fisher, run_mwu, run_roc, subset_string, temp_folder):
	start = time.time()
	if subset_string is None:
		data = read_input(input_file, mode='r+')
		X = data.X[:]
		output_file = input_file
	else:
		attr, value = subset_string.split(':')
		data = read_input(input_file, mode='a')
		data = data[data.obs[attr] == value].copy()
		X = data.X
		import os
		output_file = os.path.splitext(output_excel_file)[0] + '.h5ad'

	end = time.time()
	print("{0} is loaded. Time spent = {1:.2f}s.".format(input_file, end - start))
	__markers(data, X, labels = labels, n_jobs = n_jobs, run_fisher = run_fisher, run_mwu = run_mwu, run_roc = run_roc, temp_folder=temp_folder)
	data.write(output_file)

	print("Differential expression results are written back to h5ad file.")

	write_results_to_excel(output_excel_file, data.var, alpha=alpha)
