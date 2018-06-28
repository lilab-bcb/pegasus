#!/usr/bin/env python

import anndata
import numpy as np
import pandas as pd
from natsort import natsorted
from scipy.stats import f_oneway

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
	data = anndata.read_h5ad(input_file, backed = 'r')
	if show_attributes:
		print("Available sample attributes in input dataset: {0}".format(', '.join(data.obs.columns.values)))
	if show_gene_attributes:
		print("Available gene attributes in input dataset: {0}".format(', '.join(data.var.columns.values)))
	if not show_values_for_attributes is None:
		for attr in show_values_for_attributes.split(','):
			print("Available values for attribute {0}: {1}.".format(attr, ', '.join(np.unique(data.obs[attr]))))
