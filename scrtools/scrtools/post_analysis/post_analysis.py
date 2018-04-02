#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc

def search_gene(adata, gene, test = 'fisher', direction = 'up'):
	ncluster = adata.obs['louvain_groups'].nunique()
	results = []
	for cid in range(ncluster):
		nstr = 'de_{test}_{cid}_{direction}'.format(test = test, cid = cid, direction = direction)
		idx = np.isin(adata.uns[nstr + '_genes'], gene).nonzero()[0]
		assert idx.size <= 1
		if idx.size == 0:
			continue
		results.append(pd.DataFrame(adata.uns[nstr + '_stats'][idx[0] : idx[0]+1], index = pd.Index([str(cid + 1)], name = "Cluster ID")))
	if len(results) > 0:
		results = pd.concat(results)
	else:
		results = None
	return results

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
