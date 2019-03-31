from .Base import Base
from ..pipeline import run_pipeline

class SubClustering(Base):
    """
Run scCloud to obtain subclusters.

Usage:
  scCloud subcluster [options] --subset-selection <subset-selection>... <input_file> <output_name>
  scCloud subcluster -h

Arguments:
  input_file             Single cell data with clustering done in h5ad format.
  output_name            Output file name. All outputs will use it as the prefix.

Options:
  --subset-selection <subset-selection>...         Specify which cells will be included in the subcluster analysis. Each <subset_selection> string takes the format of 'attr:value,...,value', which means select cells with attr in the values. If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings.

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
  --diffmap-full-speed                             For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

  --calculate-pseudotime <roots>                   Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  --run-louvain                                    Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --louvain-affinity <affinity>                    Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

  --run-approximated-louvain                       Run approximated louvain clustering algorithm.
  --approx-louvain-ninit <number>                  Number of Kmeans tries. [default: 20]
  --approx-louvain-nclusters <number>              Number of clusters for Kmeans initialization. [default: 30]
  --approx-louvain-resolution <resolution>.        Resolution parameter for louvain. [default: 1.3]

  --run-tsne                                       Run multi-core t-SNE for visualization.
  --run-net-tsne                                   Run net tSNE for visualization.
  --run-fitsne                                     Run FIt-SNE for visualization.
  --run-net-fitsne                                 Run net FIt-SNE for visualization.
  --tsne-perplexity <perplexity>                   t-SNE's perplexity parameter. [default: 30]

  --run-umap                                       Run umap for visualization.
  --run-net-umap                                   Run net umap for visualization.
  --umap-K <K>                                     K neighbors for umap. [default: 15]
  --umap-min-dist <number>                         Umap parameter. [default: 0.1]
  --umap-spread <spread>                           Umap parameter. [default: 1.0]

  --run-fle                                        Run force-directed layout embedding.
  --run-net-fle                                    Run net FLE.
  --fle-K <K>                                      K neighbors for building graph for FLE. [default: 50]
  --fle-target-change-per-node <change>            Target change per node to stop forceAtlas2. If we cannot reach the target change within 5,000 iterations, we stop the algorithm as well. [default: 2.0]

  --net-knn-indices <string>                       kNN indices used for net-related visualization, can be either knn_indices or diffmap_knn_indices. [default: diffmap_knn_indices]
  --net-first-K <first_K>                          Cover <first_K> neighbors for net-related visualization. [default: 5]

  -h, --help                                       Print out help information.

Outputs:
  output_name.h5ad              Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, t-SNE and diffusion map coordinates are stored in the 'obsm' field.
  output_name.loom              Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.
  
Examples:
  scCloud subcluster -p 20 --correct-batch-effect --subset-selection louvain_labels:3,6 --subset-selection Condition:CB_nonmix --run-tsne --run-louvain manton_bm.h5ad manton_bm_subset
    """

    def execute(self):
        kwargs = {
            'processed' : True,
            'subcluster' : True,
            'cite_seq' : False,
            
            'subset_selections' : self.args['--subset-selection'],

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
            'diffmap_full_speed' : self.args['--diffmap-full-speed'],

            'run_louvain' : self.args['--run-louvain'],
            'louvain_resolution' : float(self.args['--louvain-resolution']),
            'louvain_affinity' : self.args['--louvain-affinity'],

            'run_approx_louvain' : self.args['--run-approximated-louvain'],
            'approx_louvain_ninit' : int(self.args['--approx-louvain-ninit']),
            'approx_louvain_nclusters' : int(self.args['--approx-louvain-nclusters']),
            'approx_louvain_resolution' : float(self.args['--approx-louvain-resolution']),

            'run_tsne' : self.args['--run-tsne'],
            'run_net_tsne' : self.args['--run-net-tsne'],
            'run_fitsne' : self.args['--run-fitsne'],
            'run_net_fitsne' : self.args['--run-net-fitsne'],
            'tsne_perplexity' : float(self.args['--tsne-perplexity']),

            'run_umap' : self.args['--run-umap'],
            'run_net_umap' : self.args['--run-net-umap'],
            'umap_K' : int(self.args['--umap-K']),
            'umap_min_dist' : float(self.args['--umap-min-dist']),
            'umap_spread' : float(self.args['--umap-spread']),

            'run_fle' : self.args['--run-fle'],
            'run_net_fle' : self.args['--run-net-fle'],
            'fle_K' : int(self.args['--fle-K']),
            'fle_target_change_per_node' : float(self.args['--fle-target-change-per-node']),

            'knn_indices' : self.args['--net-knn-indices'],
            'first_K' : int(self.args['--net-first-K']),

            'pseudotime' : self.split_string(self.args['--calculate-pseudotime'])
        }

        run_pipeline(self.args['<input_file>'], self.args['<output_name>'], **kwargs)
