import os
import numpy as np
import pandas as pd
import scanpy.api as sc

from . import batch_correction



def subcluster(input_h5ad_file, cluster_ids, output_name, **kwargs):
	output_name = os.path.abspath(output_name)
	output_dir = os.path.dirname(output_name)
	output_base = '.' + os.path.basename(output_name)

	sc.settings.verbosity = 3
	sc.settings.writedir = sc.settings.figdir = output_dir + '/'

	# Select a cluster to do sub-cluster
	adata = sc.read(input_h5ad_file)
	obs_index = np.isin(adata.obs['louvain_labels'], cluster_ids)
	adata._inplace_subset_obs(obs_index)

	# Select variable genes
	# ~ 10^3, for 10^5, will include majority of expressed genes
	# filter_result contains a boolean list (gene_subset) within robust genes
	filter_result = batch_correction.filter_genes_dispersion(adata, kwargs['correct_batch'])
	if kwargs['diagnosis']:
		sc.pl.filter_genes_dispersion(filter_result, save=output_base)
		print("Plotted variable gene plot.")

	adata_c = batch_correction.collect_variable_gene_matrix(adata, filter_result.gene_subset)

	keys = ['louvain_labels']
	if kwargs['key'] is not None:
		keys.append(kwargs['key'])

	# Correct batch effects
	if kwargs['correct_batch']:
		batch_correction.correct_batch_effects(adata_c)
	
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
	sc.tl.tsne(adata_c, n_jobs = kwargs['nthreads'])
	if kwargs['diagnosis']:
		sc.pl.tsne(adata_c, save=output_base)
		print("Plotted tSNE plot.")

	adata_c.write(output_name + '_var.h5ad')
	print("Written tsne results.")

	# Louvain clustering
	sc.tl.louvain(adata_c, n_neighbors = 50, resolution = kwargs['resolution'], n_jobs = kwargs['nthreads'])
	adata_c.obs['louvain_labels'] = [str(int(x) + 1) for x in adata_c.obs['louvain_groups']]

	adata_c.write(output_name + '_var.h5ad')
	if kwargs['output_loom']:
		adata_c.write_loom(output_name + '_var.loom')
	print("Written louvain cluster results.")

	# Copy clustering information to the big matrix
	adata.obs['louvain_groups'] = adata_c.obs['louvain_groups']
	adata.obs['louvain_labels'] = adata_c.obs['louvain_labels']
	adata.obsm['X_tsne'] = adata_c.obsm['X_tsne']
 
	adata.write(output_name + '.h5ad')
	if kwargs['output_loom']:
		adata.write_loom(output_name + '.loom')
	print("Written full data.")

	# Generate tSNE plots
	if kwargs['legend_on_data']:
		sc.pl.tsne(adata_c, color = keys, save = output_base, legend_loc = "on data")
	else:
		sc.pl.tsne(adata_c, color = keys, save = output_base, legend_fontsize = 10)
