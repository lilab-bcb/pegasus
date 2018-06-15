#!/usr/bin/env python

import numpy as np
import pandas as pd
from natsort import natsorted
from scipy.stats import f_oneway

def search_genes(data, gene_list, measure = 'percentage'):
	ncluster =  sum([1 for x in data.var.columns if x.startswith('mean_log_expression_')])
	columns = ["{}_{}".format(measure, x + 1) for x in range(ncluster)]
	return data.var.loc[gene_list, columns]

def search_de_genes(data, gene_list, test = 'fisher', thre = 1.5):
	ngene = len(gene_list)
	ncluster = sum([1 for x in data.var.columns if x.startswith('mean_log_expression_')])
	results = np.zeros((ngene, ncluster), dtype = np.dtype('U4'))
	columns = [str(x + 1) for x in range(ncluster)]
	df_de = data.var.loc[gene_list, [test + "_qval_" + x for x in columns]]
	if test == 'fisher':
		df_fc = data.var.loc[gene_list, ["percentage_fold_change_" + x for x in columns]]
	else:
		df_fc = np.exp(data.var.loc[gene_list, ["log_fold_change_" + x for x in columns]])
	results[:] = '?'
	results[np.isnan(df_de)] = 'NaN'
	results[(df_de <= 0.05).values & (df_fc > 1.0).values] = '+'
	results[(df_de <= 0.05).values & (df_fc >= thre).values] = '++'
	results[(df_de <= 0.05).values & (df_fc < 1.0).values] = '-'
	results[(df_de <= 0.05).values & (df_fc <= 1.0 / thre).values] = '--'
	df = pd.DataFrame(results, index = gene_list, columns = columns)
	return df
