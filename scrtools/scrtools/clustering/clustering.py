import os
import time
import numpy as np
import pandas as pd
import scanpy.api as sc 
import anndata
import tables
import xlsxwriter
from scipy.sparse import csr_matrix
from collections import Counter

from . import batch_correction, diffusion_map, dimension_reduction



def cluster(input_name, output_name, **kwargs):
	output_name = os.path.abspath(output_name)
	output_dir = os.path.dirname(output_name)
	output_base = '.' + os.path.basename(output_name)

	sc.settings.verbosity = 3
	sc.settings.writedir = sc.settings.figdir = output_dir + '/'

	if not kwargs['input_h5ad']:
		adata = 
	else:
		adata = 
	
	update_var_names(adata, kwargs['genome'])



	# Obtain normalized and logged gene expression matrix
	batch_correction.normalization(adata, kwargs['norm_count'])
	adata.X = adata.X

	# Augment key barcode annotations for batch correction
	if not kwargs['input_h5ad']:
		df = pd.read_csv(input_name + '.attr.csv', header = 0, index_col = 0, dtype = str)
		df = df.loc[adata.obs['Channel'], kwargs['attrs']]
		df.index = adata.obs.index
		adata.obs[kwargs['attrs']] = df

	if kwargs['correct_batch']:
		if len(kwargs['groupby']) == 0:
			adata.obs['Group'] = 'one_group'
		else:
			adata.obs['Group'] = adata.obs[kwargs['groupby']].apply(lambda x: '+'.join([str(y) for y in x]), axis = 1)

	# Calculate batch correction adjustment matrices
	if kwargs['correct_batch']:
		batch_correction.estimate_adjustment_matrices(adata)

	adata.write(output_name + '.h5ad')
	print("Written log normalized data.")

	# Select variable genes
	# ~ 10^3, for 10^5, will include majority of expressed genes
	# filter_result contains a boolean list (gene_subset) within robust genes
	filter_result = batch_correction.filter_genes_dispersion(adata, kwargs['correct_batch'])
	# if kwargs['diagnosis']:
	# 	sc.pl.filter_genes_dispersion(filter_result, save=output_base)
	# 	print("Plotted variable gene plot.")

	adata_c = batch_correction.collect_variable_gene_matrix(adata, filter_result.gene_subset)


	keys = ['louvain_labels']
	if kwargs['key'] is not None:
		keys.append(kwargs['key'])

	# Correct batch effects
	if kwargs['correct_batch']:
		batch_correction.correct_batch_effects(adata_c)
		from scipy.sparse import issparse
		assert issparse(adata_c.X)
		print("PASS!")

	sc.pp.scale(adata_c, max_value=10) # if x > max_value, x = max_value; since we have many zeros, negative values will have much smaller magnitudes than positive values.

	# PCA analysis
	sc.tl.pca(adata_c)
	adata_c.obsm['X_pca'] *= -1
	
	if kwargs['diagnosis']:
		sc.pl.pca_scatter(adata_c, right_margin=0.2, save=output_base)
		sc.pl.pca_variance_ratio(adata_c, save=output_base)
		sc.pl.pca_loadings(adata_c, save=output_base)
		print("Plotted PCA plots.")

	adata_c.write(output_name + '_var.h5ad')
	print("Written pca results.")

	# tSNE analysis
	# sc.tl.tsne(adata_c, n_jobs = kwargs['nthreads'])

	# W = diffusion_map.calculate_symmetric_knn_graph(adata_c.obsm['X_pca'], kwargs['nthreads'])
	# Phi_t, t, Phi = diffusion_map.calculate_diffusion_map(W, n_components = 50, sigma = 0.01)
	# print("t = {}.".format(t))
	# adata_c.obsm['X_dmap'] = Phi
	# from sklearn.decomposition import PCA
	# pca = PCA(n_components = 2, random_state = 0)
	# X_pca = pca.fit_transform(adata_c.obsm['X_dmap'])
	# print(X_pca)
	# adata_c.obsm['X_tsne'] = Phi[:, 0:2]

	# dimension_reduction.run_tsne(adata_c, 'X_dmap', kwargs['nthreads'])
	# if kwargs['diagnosis']:
	# 	sc.pl.tsne(adata_c, save=output_base)
	# 	print("Plotted tSNE plot.")

	# adata_c.write(output_name + '_var.h5ad')
	# print("Written tsne results.")

	# Louvain clustering
	sc.tl.louvain(adata_c, n_neighbors = 50, resolution = kwargs['resolution'], n_jobs = kwargs['nthreads'])
	adata_c.obs['louvain_labels'] = [str(int(x) + 1) for x in adata_c.obs['louvain_groups']]

	adata_c.write(output_name + '_var.h5ad')
	if kwargs['output_loom']:
		adata_c.write_loom(output_name + '_var.loom')
	print("Written louvain cluster results.")

	# dimension_reduction.run_tsne(adata_c, 'X_diffmap', kwargs['nthreads'])
	from sklearn.decomposition import PCA
	pca = PCA(n_components = 2, random_state = 0)
	X_pca = pca.fit_transform(adata_c.obsm['X_diffmap'])
	adata_c.obsm['X_tsne'] = X_pca

	# Copy clustering information to the big matrix
	# adata.obs['louvain_groups'] = adata_c.obs['louvain_groups']
	# adata.obs['louvain_labels'] = adata_c.obs['louvain_labels']
	# adata.obsm['X_tsne'] = adata_c.obsm['X_tsne']
 
	# adata.write(output_name + '.h5ad')
	# if kwargs['output_loom']:
	# 	adata.write_loom(output_name + '.loom')
	# print("Written full data.")

	# Generate tSNE plots
	if kwargs['legend_on_data']:
		sc.pl.tsne(adata_c, color = keys, save = output_base, legend_loc = "on data")
	else:
		sc.pl.tsne(adata_c, color = keys, save = output_base, legend_fontsize = 10)
