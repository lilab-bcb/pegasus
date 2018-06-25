#!/usr/bin/env python

import numpy as np
import pandas as pd
from natsort import natsorted

def search_genes(adata, gene_list, measure = 'percentage'):
	if not isinstance(gene_list, np.ndarray):
		gene_list = np.array(gene_list)
	gene_list = gene_list[np.isin(gene_list, adata.var.index.values)]
	ncluster = adata.obs['louvain_labels'].nunique()
	columns = natsorted(adata.obs['louvain_labels'].unique())
	gene2loc = pd.Series(range(adata.shape[1]), index = adata.var.index)
	gidx = gene2loc[gene_list].values
	results = np.zeros((gidx.size, ncluster))
	for i, cid in enumerate(columns):
		idx = (adata.obs['louvain_labels'] == cid).values
		mat = adata.X[idx][:,gidx]
		if measure == 'percentage':
			results[:, i] = mat.getnnz(axis = 0) * 100.0 / mat.shape[0]
		else:
			assert measure == 'mean_log_expression'
			results[:, i] = mat.sum(axis = 0).A1 / mat.shape[0]
	df = pd.DataFrame(results, index = gene_list, columns = columns)
	return df

def search_de_gene(adata, gene, test = 'fisher', direction = 'up'):
	ncluster = adata.obs['louvain_labels'].nunique()
	clusters = [str(x + 1) for x in range(ncluster)]
	results = []
	for cid in clusters:
		nstr = 'de_{test}_{cid}_{direction}'.format(test = test, cid = cid, direction = direction)
		idx = np.isin(adata.uns[nstr + '_genes'], gene).nonzero()[0]
		assert idx.size <= 1
		if idx.size == 0:
			continue
		results.append(pd.DataFrame(adata.uns[nstr + '_stats'][idx[0] : idx[0]+1], index = pd.Index([cid], name = "Cluster ID")))
	if len(results) > 0:
		results = pd.concat(results)
	else:
		results = None
	return results

def search_de_genes(adata, gene_list, test = 'fisher', thre = 1.5):
	if not isinstance(gene_list, np.ndarray):
		gene_list = np.array(gene_list)
	fc = 'fold_change' if test == 'fisher' else 'log_fold_change'
	ngene = gene_list.size
	ncluster = adata.obs['louvain_labels'].nunique()
	columns = [str(x + 1) for x in range(ncluster)]
	results = np.zeros((ngene, ncluster), dtype = np.dtype('U3'))
	for i, cid in enumerate(columns):
		results[:, i] = '?'
		colgene = 'de_{test}_{cid}_up_genes'.format(test = test, cid = cid)
		colstat = 'de_{test}_{cid}_up_stats'.format(test = test, cid = cid)
		idx = np.isin(gene_list, adata.uns[colgene])
		results[idx, i] = '+'
		if idx.sum() > 0:	
			idx_loc = idx.nonzero()[0]
			gene2loc = pd.Series(range(adata.uns[colgene].size), index = adata.uns[colgene])
			pos_list = gene2loc[gene_list[idx_loc]].values
			idx2 = adata.uns[colstat][pos_list][fc] >= thre
			results[idx_loc[idx2], i] = '++'
		colgene = 'de_{test}_{cid}_down_genes'.format(test = test, cid = cid)
		colstat = 'de_{test}_{cid}_down_stats'.format(test = test, cid = cid)
		idx = np.isin(gene_list, adata.uns[colgene])
		results[idx, i] = '-'
		if idx.sum() > 0:	
			idx_loc = idx.nonzero()[0]
			gene2loc = pd.Series(range(adata.uns[colgene].size), index = adata.uns[colgene])
			pos_list = gene2loc[gene_list[idx_loc]].values
			idx2 = adata.uns[colstat][pos_list][fc] <= 1.0 / thre
			results[idx_loc[idx2], i] = '--'
	df = pd.DataFrame(results, index = gene_list, columns = columns)
	return df

def calc_gene_stat(adata, clust_id):
	clust_id = str(clust_id)
	idx = np.asarray(adata.obs['louvain_labels'] == clust_id)
	mat = adata.X[idx,:]
	df = pd.DataFrame({'percentage' : mat.getnnz(axis = 0) * 100.0 / mat.shape[0],
					   'mean_log_expression' : mat.sum(axis = 0).A1 / mat.shape[0]},
					   index = pd.Index(adata.var.index, name = 'gene'),
					   columns = ['percentage', 'mean_log_expression'])
	df.sort_values(by = 'percentage', ascending = False, inplace = True)
	return df
