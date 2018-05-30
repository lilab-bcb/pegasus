from .Base import Base
from ..pipeline import run_pipeline

class Clustering(Base):
    """
Run scrtools.pipeline to obtain top-level clusters.

Usage:
  scrtools cluster [options] <input_file> <output_name>
  scrtools cluster -h

Arguments:
  input_file       Input file in either 10x format or h5ad format.
  output_name      Output name. output_name.h5ad and output_name.log will be generated.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --genome <genome>                                Genome name. [default: GRCh38]
  --processed                                      Input file is processed and thus no PCA & diffmap will be run.

  --output-filtration-results <spreadsheet>        Output filtration results into <spreadsheet>.
  --output-loom                                    Output loom-formatted file.
  --correct-batch-effect                           Correct for batch effects.
  --batch-group-by <expression>                    Group batches according to <expression>. If <expression> is None, assume all channels are of one group.
  
  --min-genes <number>                             Only keep cells with at least <number> of genes. [default: 500]
  --max-genes <number>                             Only keep cells with less than <number> of genes. [default: 6000]
  --mito-prefix <prefix>                           Prefix for mitochondrial genes. [default: MT-]
  --percent-mito <ratio>                           Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]
  --gene-percent-cells <ratio>                     Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]

  --counts-per-cell-after <number>                 Total counts per cell after normalization. [default: 1e5]
  
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
  scrtools cluster -p 20 --correct-batch-effect --run-louvain --run-tsne manton_bm_10x.h5 manton_bm
    """

    def execute(self):
        kwargs = {
            'n_jobs' : int(self.args['--threads']),
            'genome' : self.args['--genome'],

            'processed' : self.args['--processed'],
            'subcluster' : False,

            'filt_xlsx' : self.args['--output-filtration-results'],
            'output_loom' : self.args['--output-loom'],
            'batch_correction' : self.args['--correct-batch-effect'],
            'group_attribute' : self.args['--batch-group-by'],

            'min_genes' : int(self.args['--min-genes']),
            'max_genes' : int(self.args['--max-genes']),
            'mito_prefix' : self.args['--mito-prefix'],
            'percent_mito' : float(self.args['--percent-mito']),
            'percent_cells' : float(self.args['--gene-percent-cells']),

            'norm_count' : float(self.args['--counts-per-cell-after']),

            'random_state' : int(self.args['--random-state']),
            
            'pca_key' : 'X_pca' if not self.args['--run-uncentered-pca'] else 'X_rpca',
            'select_variable_genes' : not self.args['--no-variable-gene-selection'],
            'submat_to_dense' : not self.args['--no-submat-to-dense'],
            'nPC' : int(self.args['--nPC']),

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
