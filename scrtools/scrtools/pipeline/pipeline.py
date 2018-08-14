import os
from .. import tools 

def run_pipeline(input_file, output_name, **kwargs):
	# load input data
	is_raw = not kwargs['processed']
	adata = tools.read_input(input_file, genome = kwargs['genome'])

	# preprocessing
	if is_raw:
		# make gene names unique
		tools.update_var_names(adata, kwargs['genome'])
		# filter out low quality cells/genes
		tools.filter_data(adata, mito_prefix = kwargs['mito_prefix'], filt_xlsx = kwargs['filt_xlsx'], min_genes = kwargs['min_genes'], max_genes = kwargs['max_genes'], percent_mito = kwargs['percent_mito'], percent_cells = kwargs['percent_cells'])
		# normailize counts and then transform to log space
		tools.log_norm(adata, kwargs['norm_count'])
		# estimate bias factors
		if kwargs['batch_correction']:
			tools.set_group_attribute(adata, kwargs['group_attribute'])
			tools.estimate_adjustment_matrices(adata)
	elif kwargs['subcluster']:
		adata = tools.get_anndata_for_subclustering(adata, kwargs['subset_selections'])
		is_raw = True # get submat and then set is_raw to True


	if kwargs['output_seurat_compatible']:
		assert is_raw and kwargs['select_variable_genes'] and kwargs['submat_to_dense']	


	# dimension reduction --- select variable genes or not
	pca_key = kwargs['pca_key']
	if is_raw:	
		if kwargs['select_variable_genes']:
			filter_result = tools.filter_genes_dispersion(adata, kwargs['batch_correction'])
			adata_c = tools.collect_variable_gene_matrix(adata, filter_result.gene_subset)
			if kwargs['submat_to_dense']:
				adata_c.X = adata_c.X.toarray()
			if kwargs['batch_correction']:
				tools.correct_batch_effects(adata_c)
		
			# dimension reduction
			if pca_key == 'X_pca':
				tools.run_pca(adata_c, nPC = kwargs['nPC'], random_state = kwargs['random_state'])
			else:
				tools.run_rpca(adata_c, nPC = kwargs['nPC'], random_state = kwargs['random_state'])
			adata.obsm[pca_key] = adata_c.obsm[pca_key]
		else:
			assert pca_key == 'X_rpca'
			if kwargs['batch_correction']:
				tools.correct_batch_effects(adata)
			tools.run_rpca(adata, nPC = kwargs['nPC'], random_state = kwargs['random_state'])
	else:
		assert pca_key in adata.obsm.keys()

	# diffusion map
	if is_raw:
		tools.run_diffmap(adata, pca_key, n_jobs = kwargs['n_jobs'], n_components = kwargs['nDC'], alpha = kwargs['diffmap_alpha'], K = kwargs['diffmap_K'], random_state = kwargs['random_state'], full_speed = kwargs['diffmap_full_speed'])
	else:
		assert 'X_diffmap' in adata.obsm.keys()

	# clustering
	if kwargs['run_approx_louvain']:
		tools.run_approximated_louvain(adata, 'X_diffmap', n_jobs = kwargs['n_jobs'], resolution = kwargs['approx_louvain_resolution'], random_state = kwargs['random_state'], n_clusters = kwargs['approx_louvain_nclusters'], n_init = kwargs['approx_louvain_ninit'])
	if kwargs['run_kmeans']:
		tools.run_kmeans(adata, 'X_diffmap', kwargs['kmeans_n_clusters'], n_jobs = kwargs['n_jobs'], random_state = kwargs['random_state'])
	if kwargs['run_louvain']:
		tools.run_louvain(adata, affinity = kwargs['louvain_affinity'], resolution = kwargs['louvain_resolution'], random_state = kwargs['random_state'])
	if kwargs['run_hdbscan']:
		tools.run_hdbscan(adata, 'X_diffmap', n_jobs = kwargs['n_jobs'], min_cluster_size = kwargs['hdbscan_min_cluster_size'], min_samples = kwargs['hdbscan_min_samples'])

	# visualization
	if kwargs['run_tsne']:
		tools.run_tsne(adata, pca_key, n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'])
	if kwargs['run_fitsne']:
		tools.run_fitsne(adata, pca_key, n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'])
	if kwargs['run_umap']:
		if kwargs['run_umap_on_diffmap']:
			tools.run_umap(adata, 'X_diffmap', n_neighbors = kwargs['umap_K'], min_dist = kwargs['umap_min_dist'], spread = kwargs['umap_spread'], random_state = kwargs['random_state'])
			adata.obsm['X_umap_diffmap'] = adata.obsm['X_umap']
		tools.run_umap(adata, pca_key, n_neighbors = kwargs['umap_K'], min_dist = kwargs['umap_min_dist'], spread = kwargs['umap_spread'], random_state = kwargs['random_state'])
	if kwargs['run_fle']:
		tools.run_force_directed_layout(adata, output_name, affinity = kwargs['fle_affinity'], n_jobs = kwargs['n_jobs'], K = kwargs['fle_K'], n_steps = kwargs['fle_n_steps'])	
	
	# calculate diffusion-based pseudotime from roots
	if kwargs['pseudotime'] is not None:
		assert 'X_diffmap' in adata.obsm.keys()
		tools.run_pseudotime_calculation(adata, kwargs['pseudotime'])

	adata.write(output_name + ".h5ad")

	if kwargs['output_seurat_compatible']:
		adata_c.obs = adata.obs
		adata_c.obsm = adata.obsm
		adata_c.uns = adata.uns
		adata_c.raw = adata
		adata_c.write(output_name + ".seurat.h5ad")

	if kwargs['output_loom']:
		adata.write_loom(output_name + ".loom")

	print("Results are written.")
