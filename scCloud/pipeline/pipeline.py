import numpy as np
import anndata
from scipy.sparse import csr_matrix, hstack
from scCloud import tools, cite_seq

def run_pipeline(input_file, output_name, **kwargs):
	is_raw = not kwargs['processed']

	if 'seurat_compatible' not in kwargs:
		kwargs['seurat_compatible'] = False

	# load input data
	if not kwargs['cite_seq']:
		adata = tools.read_input(input_file, genome = kwargs['genome'], mode = ('a' if (is_raw or kwargs['subcluster']) else 'r+'), select_singlets = kwargs['select_singlets'])
	else:
		data_dict = tools.read_input(input_file, genome = kwargs['genome'], return_a_dict = True, select_singlets = kwargs['select_singlets'])
		assert len(data_dict) == 2
		adata = cdata = None
		for genome, data in data_dict.items():
			if genome.startswith('CITE_Seq'):
				cdata = data
			else:
				adata = data
		assert adata is not None and cdata is not None
	print("Inputs are loaded.")

	if kwargs['seurat_compatible']:
		assert is_raw and kwargs['select_variable_genes'] and kwargs['submat_to_dense']

	# preprocessing
	if is_raw:
		# make gene names unique
		tools.update_var_names(adata)
		# filter out low quality cells/genes
		tools.filter_data(adata, output_filt = kwargs['output_filt'], plot_filt = kwargs['plot_filt'], plot_filt_figsize = kwargs['plot_filt_figsize'], \
			mito_prefix = kwargs['mito_prefix'], min_genes = kwargs['min_genes'], max_genes = kwargs['max_genes'], min_umis = kwargs['min_umis'], max_umis = kwargs['max_umis'], \
			percent_mito = kwargs['percent_mito'], percent_cells = kwargs['percent_cells'], min_genes_on_raw = kwargs['min_genes_on_raw'])
		if kwargs['seurat_compatible']:
			raw_data = adata.copy() # raw as count
		# normailize counts and then transform to log space
		tools.log_norm(adata, kwargs['norm_count'])
		# estimate bias factors
		if kwargs['batch_correction']:
			tools.set_group_attribute(adata, kwargs['group_attribute'])
			tools.estimate_adjustment_matrices(adata)
	elif kwargs['subcluster']:
		adata = tools.get_anndata_for_subclustering(adata, kwargs['subset_selections'])
		is_raw = True # get submat and then set is_raw to True

	# dimension reduction --- select variable genes or not
	if is_raw:
		if kwargs['select_hvg']:
			tools.select_highly_variable_genes(adata, kwargs['batch_correction'], flavor = kwargs['hvg_flavor'], n_top = kwargs['hvg_ngenes'], plot_hvg_fig = kwargs['plot_hvg'], n_jobs = kwargs['n_jobs'], benchmark_time = kwargs.get('benchmark_time', False))
		adata_c = tools.collect_highly_variable_gene_matrix(adata) # select hvg matrix and convert to dense
		if kwargs['batch_correction']:
			tools.correct_batch_effects(adata_c)
		tools.run_pca(adata_c, nPC = kwargs['nPC'], random_state = kwargs['random_state'])
		adata.obsm['X_pca'] = adata_c.obsm['X_pca']
	else:
		assert 'X_pca' in adata.obsm.keys()


	# diffusion map
	if is_raw:
		tools.run_diffmap(adata, 'X_pca', n_jobs = kwargs['n_jobs'], n_components = kwargs['nDC'], alpha = kwargs['diffmap_alpha'], K = kwargs['diffmap_K'], solver = kwargs['diffmap_solver'], random_state = kwargs['random_state'], full_speed = kwargs['diffmap_full_speed'])

		import time
		start_time = time.time()
		tools.get_kNN(adata, 'X_diffmap', kwargs['diffmap_K'], n_jobs = kwargs['n_jobs'], random_state = kwargs['random_state'], full_speed = kwargs['diffmap_full_speed'])
		end_time = time.time()
		print("KNN for diffusion components is finished. Time spent = {:.2f}s.".format(end_time - start_time))

		adata.obsm['X_diffmap_pca'] = tools.reduce_diffmap_to_3d(adata.obsm['X_diffmap'], random_state = kwargs['random_state'])
	else:
		assert 'X_diffmap' in adata.obsm.keys()

	# calculate kBET
	if ('kBET' in kwargs) and kwargs['kBET']:
		stat_mean, pvalue_mean, accept_rate = tools.calc_kBET(adata, kwargs['kBET_batch'], K = kwargs['kBET_K'], alpha = kwargs['kBET_alpha'], n_jobs = kwargs['n_jobs'])
		print("kBET stat_mean = {:.2f}, pvalue_mean = {:.4f}, accept_rate = {:.2%}.".format(stat_mean, pvalue_mean, accept_rate))
		# if kwargs['kBJSD']:
		# 	print("kBJSD mean = {:.4f}".format(tools.calc_kBJSD(adata, kwargs['kBET_batch'], K = kwargs['kBET_K'], n_jobs = kwargs['n_jobs'])))

	# clustering
	if kwargs['run_approx_louvain']:
		tools.run_approximated_louvain(adata, 'X_' + kwargs['approx_louvain_basis'], affinity = kwargs['approx_louvain_affinity'], resolution = kwargs['approx_louvain_resolution'], n_clusters = kwargs['approx_louvain_nclusters'], \
			n_init = kwargs['approx_louvain_ninit'], n_jobs = kwargs['n_jobs'], random_state = kwargs['random_state'], temp_folder = kwargs['temp_folder'], class_label = 'approx_louvain_labels')

	if kwargs['run_approx_leiden']:
		tools.run_approximated_leiden(adata, 'X_' + kwargs['approx_leiden_basis'], affinity = kwargs['approx_leiden_affinity'], resolution = kwargs['approx_leiden_resolution'], n_clusters = kwargs['approx_leiden_nclusters'], \
			n_init = kwargs['approx_leiden_ninit'], n_jobs = kwargs['n_jobs'], random_state = kwargs['random_state'], temp_folder = kwargs['temp_folder'], class_label = 'approx_leiden_labels')

	if kwargs['run_louvain']:
		tools.run_louvain(adata, affinity = kwargs['louvain_affinity'], resolution = kwargs['louvain_resolution'], random_state = kwargs['random_state'], class_label = kwargs['louvain_class_label'])

	if kwargs['run_leiden']:
		tools.run_leiden(adata, affinity = kwargs['leiden_affinity'], resolution = kwargs['leiden_resolution'], n_iter = kwargs['leiden_niter'], random_state = kwargs['random_state'], class_label = kwargs['leiden_class_label'])



	# visualization
	if kwargs['run_net_tsne']:
		selected = tools.select_cells(adata.uns['knn_distances'], kwargs['net_ds_frac'], K = kwargs['net_ds_K'], alpha = kwargs['net_ds_alpha'], random_state = kwargs['random_state'])
		tools.run_net_tsne(adata, 'X_pca', selected, n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'], net_alpha = kwargs['net_l2'], \
		                   polish_learning_frac = kwargs['net_tsne_polish_learing_frac'], polish_n_iter = kwargs['net_tsne_polish_niter'], out_basis = kwargs['net_tsne_basis'])

	if kwargs['run_net_fitsne']:
		selected = tools.select_cells(adata.uns['knn_distances'], kwargs['net_ds_frac'], K = kwargs['net_ds_K'], alpha = kwargs['net_ds_alpha'], random_state = kwargs['random_state'])
		tools.run_net_fitsne(adata, 'X_pca', selected, n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'], net_alpha = kwargs['net_l2'], \
		                   polish_learning_frac = kwargs['net_fitsne_polish_learing_frac'], polish_n_iter = kwargs['net_fitsne_polish_niter'], out_basis = kwargs['net_fitsne_basis'])

	if kwargs['run_net_umap']:
		selected = tools.select_cells(adata.uns['knn_distances'], kwargs['net_ds_frac'], K = kwargs['net_ds_K'], alpha = kwargs['net_ds_alpha'], random_state = kwargs['random_state'])
		tools.run_net_umap(adata, 'X_pca', selected, n_jobs = kwargs['n_jobs'], n_neighbors = kwargs['umap_K'], min_dist = kwargs['umap_min_dist'], spread = kwargs['umap_spread'], random_state = kwargs['random_state'], net_alpha = kwargs['net_l2'], \
		                   ds_full_speed = kwargs['net_ds_full_speed'], polish_learning_rate = kwargs['net_umap_polish_learing_rate'], polish_n_epochs = kwargs['net_umap_polish_nepochs'], out_basis = kwargs['net_umap_basis'])

	if kwargs['run_net_fle']:
		selected = tools.select_cells(adata.uns['diffmap_knn_distances'], kwargs['net_ds_frac'], K = kwargs['net_ds_K'], alpha = kwargs['net_ds_alpha'], random_state = kwargs['random_state'])
		tools.run_net_fle(adata, selected, output_name, n_jobs = kwargs['n_jobs'], K = kwargs['fle_K'], target_change_per_node = kwargs['fle_target_change_per_node'], \
			target_steps = kwargs['fle_target_steps'], is3d = kwargs['fle_3D'], random_state = kwargs['random_state'], ds_full_speed = kwargs['net_ds_full_speed'], \
			net_alpha = kwargs['net_l2'], polish_target_steps = kwargs['net_fle_polish_target_steps'], out_basis = kwargs['net_fle_basis'])


	if kwargs['run_tsne']:
		tools.run_tsne(adata, 'X_pca', n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'])

	if kwargs['run_fitsne']:
		tools.run_fitsne(adata, 'X_pca', n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'])

	if kwargs['run_umap']:
		tools.run_umap(adata, 'X_pca', n_neighbors = kwargs['umap_K'], min_dist = kwargs['umap_min_dist'], spread = kwargs['umap_spread'], random_state = kwargs['random_state'])

	if kwargs['run_fle']:
		tools.run_force_directed_layout(adata, output_name, n_jobs = kwargs['n_jobs'], K = kwargs['fle_K'], target_change_per_node = kwargs['fle_target_change_per_node'], target_steps = kwargs['fle_target_steps'], is3d = kwargs['fle_3D'], random_state = kwargs['random_state'])


	# calculate diffusion-based pseudotime from roots
	if kwargs['pseudotime'] is not None:
		assert 'X_diffmap' in adata.obsm.keys()
		tools.run_pseudotime_calculation(adata, kwargs['pseudotime'])

	# merge cite-seq data and run t-SNE
	if kwargs['cite_seq']:
		adt_matrix = np.zeros((adata.shape[0], cdata.shape[1]), dtype = 'float32')
		idx = adata.obs_names.isin(cdata.obs_names)
		adt_matrix[idx, :] = cdata[adata.obs_names[idx],].X.toarray()
		if abs(100.0 - kwargs['cite_seq_capping']) > 1e-4:
			cite_seq.capping(adt_matrix, kwargs['cite_seq_capping'])

		var_names = np.concatenate([adata.var_names, ['AD-' + x for x in cdata.var_names]])

		new_data = anndata.AnnData(X = hstack([adata.X, csr_matrix(adt_matrix)], format = 'csr'),
			obs = adata.obs,
			obsm = adata.obsm,
			uns = adata.uns,
			var = {'var_names' : var_names,
				   'gene_ids' : var_names,
				   'n_cells' : np.concatenate([adata.var['n_cells'].values, [0] * cdata.shape[1]]),
				   'percent_cells' : np.concatenate([adata.var['percent_cells'].values, [0.0] * cdata.shape[1]]),
				   'robust' : np.concatenate([adata.var['robust'].values, [False] * cdata.shape[1]]),
				   'highly_variable_genes' : np.concatenate([adata.var['highly_variable_genes'].values, [False] * cdata.shape[1]])
				  })
		new_data.obsm['CITE-Seq'] = adt_matrix
		adata = new_data
		print("ADT count matrix is attached.")

		tools.run_fitsne(adata, 'CITE-Seq', n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'], out_basis = 'citeseq_fitsne')
		print("Antibody embedding is done.")


	# write out results
	tools.write_output(adata, output_name)

	if kwargs['seurat_compatible']:
		seurat_data = adata.copy()
		seurat_data.raw = raw_data
		seurat_data.uns['scale.data'] = adata_c.X
		seurat_data.uns['scale.data.rownames'] = adata_c.var_names.values
		seurat_data.write(output_name + ".seurat.h5ad")

	if kwargs['output_loom']:
		adata.write_loom(output_name + ".loom")

	print("Results are written.")
