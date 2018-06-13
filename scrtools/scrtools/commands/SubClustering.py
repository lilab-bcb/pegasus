from .Base import Base
from ..pipeline import run_pipeline

class SubClustering(Base):
    """
Run scrtools to obtain subclusters.

Usage:
  scrtools subcluster [options] <input_file> <cluster_ids> <output_name>
  scrtools subcluster -h

Arguments:
  input_file             Single cell data with clustering done in h5ad format.
  cluster_ids            Tell us which clusters (separated by comma) you want to perform subcluster task. IDs start from 1 and must be clusters in input_file.
  output_name            Output name. output_name.h5ad and output_name.log will be generated.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --cluster-labels <labels>                        Which cluster results to use. [default: louvain_labels]
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
  --louvain-affinity <affinity>                    Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

  --run-kmeans                                     Run KMeans clustering algorithm on diffusion components.
  --kmeans-n-clusters <number>                     Target at <number> clusters for K means. [default: 20]

  --run-approximated-louvain                       Run approximated louvain clustering algorithm.
  --approx-louvain-ninit <number>                  Number of Kmeans tries. [default: 20]
  --approx-louvain-nclusters <number>              Number of clusters for Kmeans initialization. [default: 30]
  --approx-louvain-resolution <resolution>.        Resolution parameter for louvain. [default: 1.3]

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
  --fle-affinity <affinity>                        Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

  -h, --help                                       Print out help information.

Examples:
  scrtools subcluster -p 20 --correct-batch-effect manton_bm.h5ad 1 manton_bm_1
    """

    def execute(self):
        kwargs = {
            'processed' : True,
            'subcluster' : True,
            
            'cluster_labels' : self.args['--cluster-labels'],
            'cluster_ids' : self.split_string(self.args['<cluster_ids>']),

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
            'louvain_affinity' : self.args['--louvain-affinity'],

            'run_kmeans' : self.args['--run-kmeans'],
            'kmeans_n_clusters' : int(self.args['--kmeans-n-clusters']),

            'run_hdbscan' : self.args['--run-hdbscan'],
            'hdbscan_min_cluster_size' : int(self.args['--hdbscan-min-cluster-size']),
            'hdbscan_min_samples' : int(self.args['--hdbscan-min-samples']),

            'run_approx_louvain' : self.args['--run-approximated-louvain'],
            'approx_louvain_ninit' : int(self.args['--approx-louvain-ninit']),
            'approx_louvain_nclusters' : int(self.args['--approx-louvain-nclusters']),
            'approx_louvain_resolution' : float(self.args['--approx-louvain-resolution']),

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
            'fle_affinity' : self.args['--fle-affinity'],

            'pseudotime' : self.split_string(self.args['--calculate-pseudotime'])
        }

        run_pipeline(self.args['<input_file>'], self.args['<output_name>'], **kwargs)
