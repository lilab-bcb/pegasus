import os
from .. import tools 

def run_cite_seq_pipeline(input_file, output_name, **kwargs):
	# load input data
	adata = tools.read_input(input_file, genome = kwargs['genome'])
	# filter out low quality cells/genes
	tools.filter_cells_cite_seq(adata, kwargs['max_cells'])
	# normailize counts and then transform to log space
	tools.log_norm(adata, kwargs['norm_count'])

	pca_key = 'X_pca'
	adata.obsm[pca_key] = adata.X.toarray()

	# diffusion map
	tools.run_diffmap(adata, pca_key, n_jobs = kwargs['n_jobs'], n_components = kwargs['nDC'], alpha = kwargs['diffmap_alpha'], K = kwargs['diffmap_K'], random_state = kwargs['random_state'], full_speed = kwargs['diffmap_full_speed'])

	# clustering
	if kwargs['run_louvain']:
		tools.run_louvain(adata, affinity = kwargs['louvain_affinity'], resolution = kwargs['louvain_resolution'], random_state = kwargs['random_state'])
	if kwargs['run_approx_louvain']:
		tools.run_approximated_louvain(adata, 'X_diffmap', n_jobs = kwargs['n_jobs'], resolution = kwargs['approx_louvain_resolution'], random_state = kwargs['random_state'], n_clusters = kwargs['approx_louvain_nclusters'], n_init = kwargs['approx_louvain_ninit'])

	# visualization
	if kwargs['run_tsne']:
		tools.run_tsne(adata, pca_key, n_jobs = kwargs['n_jobs'], perplexity = kwargs['tsne_perplexity'], random_state = kwargs['random_state'])
		
	adata.write(output_name + ".h5ad")

	print("Results are written.")
