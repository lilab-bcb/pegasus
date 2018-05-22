import numpy as np
from scipy.sparse import issparse
import anndata
import fisher
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdr
import warnings
import pandas as pd
import scipy.stats as ss
import xlsxwriter
import threading
import time


# assume louvain_labels always start from 1 and continuous
def collect_contingency_table(data, labels = 'louvain_labels'):	
	clusts = data.obs[labels].value_counts()
	ct = np.zeros((data.var_names.size, clusts.size, 2), dtype = np.uint)
	for label, count in clusts.iteritems():
		i = int(label) - 1
		mask = np.isin(data.obs[labels], label)
		ct[:, i, 0] = data.X[mask].getnnz(axis = 0)
		ct[:, i, 1] = count - ct[:, i, 0]
	data.uns["contingency_table"] = ct

def calc_fisher_per_thread(thread_no, results, n_jobs, data, clusts, ct, total, dtypes, token):
	nclust = len(clusts)
	results[thread_no] = []
	
	for i in range(thread_no, nclust, n_jobs):
		cpt = total - ct[:, i, :]
		pvals = fisher.pvalue_npy(ct[:, i, 0], ct[:, i, 1], cpt[:, 0], cpt[:, 1])[2]
		percents = ct[:, i, 0] / (ct[:, i, 0] + ct[:, i, 1]) * 100.0

		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			fold_change = ct[:, i, 0] / (ct[:, i, 0] + ct[:, i, 1]) / (cpt[:, 0] / (cpt[:, 0] + cpt[:, 1]))
			fold_change[np.isnan(fold_change)] = 0.0

		passed, qvals = fdr(pvals, alpha = 0.05)

		grt_idx = np.logical_and(passed, fold_change > 1.0)
		results[thread_no].append(("de_fisher{0}_{1}_up_genes".format(token, clusts[i]), data.var_names[grt_idx]))
		values = list(zip(percents[grt_idx], fold_change[grt_idx], pvals[grt_idx], qvals[grt_idx]))
		results[thread_no].append(("de_fisher{0}_{1}_up_stats".format(token, clusts[i]), np.array(values, dtype = dtypes)))

		lsr_idx = np.logical_and(passed, fold_change < 1.0)
		results[thread_no].append(("de_fisher{0}_{1}_down_genes".format(token, clusts[i]), data.var_names[lsr_idx]))
		values = list(zip(percents[lsr_idx], fold_change[lsr_idx], pvals[lsr_idx], qvals[lsr_idx]))
		results[thread_no].append(("de_fisher{0}_{1}_down_stats".format(token, clusts[i]), np.array(values, dtype = dtypes)))

		print("Cluster {0} is processed.".format(clusts[i]))


# clusts is a list of strings
def fisher_test(data, clusts = None, labels = 'louvain_labels', token = "", n_jobs = 1):
	start = time.time()
	
	if "contingency_table" not in data.uns:
		collect_contingency_table(data, labels = labels)
		print("Contingency table is collected.")

	ct = data.uns["contingency_table"]

	if clusts is not None:
		ct = ct[:, [int(x) - 1 for x in clusts], :]
	else:
		clusts = [str(x + 1) for x in range(data.obs[labels].nunique())]

	total = ct.sum(axis = 1)
	dtypes = [('percentage', np.dtype('float64')), ('fold_change', np.dtype('float64')), ('pval', np.dtype('float64')), ('qval', np.dtype('float64'))]
	
	threads = [None] * n_jobs
	results = [None] * n_jobs
	for i in range(n_jobs):
		t = threading.Thread(target=calc_fisher_per_thread, args=(i, results, n_jobs, data, clusts, ct, total, dtypes, token))
		threads.append(t)
		t.start()

	for i in range(n_jobs):
		t.join()
		for key, value in results[i]:
			data.uns[key] = value

	end = time.time()
	print("Fisher's exact test is done. Time spent = {:.2f}s.".format(end - start))



def calc_t_per_thread(thread_no, results, n_jobs, data, clusts, n, v, sm1, sm2, dtypes, token):
	nclust = len(clusts)
	results[thread_no] = []
	tvar = ss.t(v)

	for i in range(thread_no, nclust, n_jobs):
		mask = np.isin(data.obs[labels], clusts[i])
		mat = data.X[mask]
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
		results[thread_no].append(("de_t{0}_{1}_up_genes".format(token, clusts[i]), data.var_names[grt_idx]))
		values = list(zip(mean1[grt_idx], log_fold_change[grt_idx], pvals[grt_idx], qvals[grt_idx]))
		results[thread_no].append(("de_t{0}_{1}_up_stats".format(token, clusts[i]), np.array(values, dtype = dtypes)))

		lsr_idx = np.logical_and(passed, log_fold_change < 1.0)
		results[thread_no].append(("de_t{0}_{1}_down_genes".format(token, clusts[i]), data.var_names[lsr_idx]))
		values = list(zip(mean1[lsr_idx], log_fold_change[lsr_idx], pvals[lsr_idx], qvals[lsr_idx]))
		results[thread_no].append(("de_t{0}_{1}_down_stats".format(token, clusts[i]), np.array(values, dtype = dtypes)))

		print("Cluster {0} is processed.".format(clusts[i]))

def t_test(data, clusts = None, labels = 'louvain_labels', token = "", n_jobs = 1):
	start = time.time()

	dtypes = [('mean_log_expression', np.dtype('float64')), ('log_fold_change', np.dtype('float64')), ('pval', np.dtype('float64')), ('qval', np.dtype('float64'))]

	if clusts is None:
		clusts = [str(x + 1) for x in range(data.obs[labels].nunique())]
	
	mask = np.isin(data.obs[labels], clusts)
	mat = data.X[mask]
	n = mat.shape[0]
	v = n - 2 # degree of freedom
	sm1 = mat.sum(axis = 0).A1 # sum of moment 1
	mat.data **= 2
	sm2 = mat.sum(axis = 0).A1 # sum of moment 2

	threads = [None] * n_jobs
	results = [None] * n_jobs
	for i in range(n_jobs):
		t = threading.Thread(target=calc_t_per_thread, args=(i, results, n_jobs, data, clusts, n, v, sm1, sm2, dtypes, token))
		threads.append(t)
		t.start()

	for i in range(n_jobs):
		t.join()
		for key, value in results[i]:
			data.uns[key] = value

	end = time.time()
	print("T test is done. Time spent = {:.2f}s.".format(end - start))



def write_results_to_excel(output_file, data, test, clusts = None, labels = 'louvain_labels', token = "", threshold = 1.5):
	if clusts is None:
		clusts = [str(x + 1) for x in range(data.obs[labels].nunique())]

	thre_kw = 'fold_change' if test == 'fisher' else 'log_fold_change'
	sort_kw = 'percentage' if test == 'fisher' else 'mean_log_expression'

	writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
	for clust_id in clusts:
		df = pd.DataFrame(data.uns["de_{0}{1}_{2}_up_stats".format(test, token, clust_id)], 
			index = pd.Index(data.uns["de_{0}{1}_{2}_up_genes".format(test, token, clust_id)], name = "gene"))
		df = df.loc[df[thre_kw] >= threshold]
		df.sort_values(by = sort_kw, ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Cluster {0}, up-regulated".format(clust_id))

		df = pd.DataFrame(data.uns["de_{0}{1}_{2}_down_stats".format(test, token, clust_id)],
			index = pd.Index(data.uns["de_{0}{1}_{2}_down_genes".format(test, token, clust_id)], name = "gene"))
		df = df.loc[df[thre_kw] <= 1.0 / threshold]
		df.sort_values(by = sort_kw, ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Cluster {0}, down-regulated".format(clust_id))
	writer.save()



def run_de_analysis(input_file, output_name, threshold, labels, n_jobs):
	data = anndata.read_h5ad(input_file)

	print("Begin t_test.")
	t_test(data, labels = labels, n_jobs = n_jobs)
	excel_file = output_name + "_de_analysis_t.xlsx"
	write_results_to_excel(excel_file, data, "t", threshold = threshold, labels = labels)
	print(excel_file + " is written.")

	print("Begin fisher_test.")
	fisher_test(data, labels = labels, n_jobs = n_jobs)
	excel_file = output_name + "_de_analysis_fisher.xlsx"
	write_results_to_excel(excel_file, data, "fisher", threshold = threshold, labels = labels)
	print(excel_file + " is written.")

	data.write(output_name + "_de.h5ad")
	print("Results are written.")
