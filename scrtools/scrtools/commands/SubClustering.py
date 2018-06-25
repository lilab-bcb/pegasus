from .Base import Base
from ..pipeline import run_pipeline

class SubClustering(Base):
    """
Run scrtools to obtain subclusters.

Usage:
  scrtools subcluster [options] [--subset-selection <subset-selection>...] <input_file> <output_name>
  scrtools subcluster -h

Arguments:
  input_file             Single cell data with clustering done in h5ad format.
  output_name            Output name. output_name.h5ad and output_name.log will be generated.

Options:
  --subset-selection <subset-selection>...         Attributes and values to perform subcluster task on. Multiple <subset_selection> strings can be specified, each <subset_selection> takes the format of 'attr:value,value'. 
  --show-attributes                                Show the available attributes in the input dataset.
  --show-values-for-attributes <attributes>        Show the available values for specified attributes in the input dataset. <attributes> should be a comma-separated list of attributes.

  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --correct-batch-effect                           Correct for batch effects for subclustering task.
  --output-loom                                    Output loom-formatted file.

  --random-state <seed>                            Random number generator seed. [default: 0]

  --run-uncentered-pca                             Run uncentered PCA.
  --no-variable-gene-selection                     Do not select variable genes.
  --no-submat-to-dense                             Do not convert variable-gene-selected submatrix to a dense matrix.
  --nPC <number>                                   Number of PCs. [default: 50]

  --nDC <number>                                   Number of diffusion components. [default: 50]
  --diffmap-alpha <alpha>                          Power parameter for diffusion-based pseudotime. [default: 0.5]
  --diffmap-K <K>                                  Number of neighbors used for constructing affinity matrix. [default: 100]

  --calculate-pseudotime <roots>                   Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a list of cell barcodes.

  --run-louvain                                    Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]

  --run-kmeans                                     Run KMeans clustering algorithm on diffusion components.
  --kmeans-n-clusters <number>                     Target at <number> clusters for K means. [default: 20]

  --run-hdbscan                                    Run hdbscan clustering algorithm on diffusion components.
  --hdbscan-min-cluster-size <number>              Minimum cluster size for hdbscan. [default: 50]
  --hdbscan-min-samples <number>                   Minimum number of samples for hdbscan. [default: 50]

  --run-tsne                                       Run multi-core tSNE for visualization.
  --tsne-perplexity <perplexity>                   tSNE's perplexity parameter. [default: 30]
  --run-fitsne                                     Run FItSNE for visualization.

  --run-umap                                       Run umap for visualization.
  --umap-on-diffmap                                Run umap on diffusion components.
  --umap-K <K>                                     K neighbors for umap. [default: 15]
  --umap-min-dist <number>                         Umap parameter. [default: 0.1]
  --umap-spread <spread>                           Umap parameter. [default: 1.0]

  --run-fle                                        Run force-directed layout embedding.
  --fle-K <K>                                      K neighbors for building graph for FLE. [default: 50]
  --fle-n-steps <nstep>                            Number of iterations for FLE. [default: 10000]

  -h, --help                                       Print out help information.

Examples:
  scrtools subcluster -p 20 --correct-batch-effect --subset_selection louvain_labels:3,6  --subset_selection Condition:CB_nonmix  manton_bm.h5ad manton_bm_subset
  scrtools subcluster -p 20 --correct-batch-effect --show_attributes manton_bm.h5ad manton_bm_tmp
  scrtools subcluster -p 20 --correct-batch-effect --show_attributes --show_values_for_attributes louvain_labels,Condition manton_bm.h5ad manton_bm_tmp
  scrtools subcluster -p 20 --correct-batch-effect --show_attributes --show_values_for_attributes louvain_labels,Condition --subset_selection Condition:CB_nonmix  manton_bm.h5ad manton_bm_subset
    """

    def execute(self):
        kwargs = {
            'processed' : True,
            'subcluster' : True,
            
            'subset_selections' : self.args['--subset-selection'],
            'show_attributes' : self.args['--show-attributes'],
            'show_values_for_attributes' : self.args['--show-values-for-attributes'],

            'n_jobs' : int(self.args['--threads']),
            'genome' : None,
            'batch_correction' : self.args['--correct-batch-effect'],
            'output_loom' : self.args['--output-loom'],

            'random_state' : int(self.args['--random-state']),
            
            'run_dimension_reduction' : True,
            'pca_key' : 'X_pca' if not self.args['--run-uncentered-pca'] else 'X_rpca',
            'select_variable_genes' : not self.args['--no-variable-gene-selection'],
            'submat_to_dense' : not self.args['--no-submat-to-dense'],
            'nPC' : int(self.args['--nPC']),

            'run_diffmap' : True,
            'nDC' : int(self.args['--nDC']),
            'diffmap_alpha' : float(self.args['--diffmap-alpha']),
            'diffmap_K' : int(self.args['--diffmap-K']),

            'run_louvain' : self.args['--run-louvain'],
            'louvain_resolution' : float(self.args['--louvain-resolution']),

            'run_kmeans' : self.args['--run-kmeans'],
            'kmeans_n_clusters' : int(self.args['--kmeans-n-clusters']),

            'run_hdbscan' : self.args['--run-hdbscan'],
            'hdbscan_min_cluster_size' : int(self.args['--hdbscan-min-cluster-size']),
            'hdbscan_min_samples' : int(self.args['--hdbscan-min-samples']),

            'run_tsne' : self.args['--run-tsne'],
            'run_fitsne' : self.args['--run-fitsne'],
            'tsne_perplexity' : float(self.args['--tsne-perplexity']),

            'run_umap' : self.args['--run-umap'],
            'run_umap_on_diffmap' : self.args['--umap-on-diffmap'],
            'umap_K' : int(self.args['--umap-K']),
            'umap_min_dist' : float(self.args['--umap-min-dist']),
            'umap_spread' : float(self.args['--umap-spread']),

            'run_fle' : self.args['--run-fle'],
            'fle_K' : int(self.args['--fle-K']),
            'fle_n_steps' : int(self.args['--fle-n-steps']),

            'pseudotime' : self.split_string(self.args['--calculate-pseudotime'])
        }

        run_pipeline(self.args['<input_file>'], self.args['<output_name>'], **kwargs)
