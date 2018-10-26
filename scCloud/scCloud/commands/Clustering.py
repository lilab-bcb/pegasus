from .Base import Base
from ..pipeline import run_pipeline

class Clustering(Base):
    """
Run scCloud.pipeline to obtain top-level clusters.

Usage:
  scCloud cluster [options] <input_file> <output_name>
  scCloud cluster -h

Arguments:
  input_file       Input file in 10x format. If first-pass analysis has been performed, but you want to run some additional analysis, you could also pass a h5ad-formatted file.
  output_name      Output file name. All outputs will use it as the prefix.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --processed                                      Input file is processed and thus no PCA & diffmap will be run.

  --genome <genome>                                A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
  
  --cite-seq                                       Data are CITE-Seq data. scCloud will perform analyses on RNA count matrix first. Then it will attach the ADT matrix to the RNA matrix with all antibody names changing to 'AD-' + antibody_name. Lastly, it will embed the antibody expression using t-SNE (the basis used for plotting is 'citeseq_tsne').

  --output-filtration-results <spreadsheet>        Output filtration results into <spreadsheet>.
  --output-seurat-compatible                       Output seurat-compatible h5ad file.
  --output-loom                                    Output loom-formatted file.
  --correct-batch-effect                           Correct for batch effects.
  --batch-group-by <expression>                    Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+...+attrn', or 'attr=value11,...,value1n_1;value21,...,value2n_2;...;valuem1,...,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',...,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,...,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.

  
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
  --diffmap-full-speed                             For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

  --calculate-pseudotime <roots>                   Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  --run-louvain                                    Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --louvain-affinity <affinity>                    Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

  --run-kmeans                                     Run KMeans clustering algorithm on diffusion components.
  --kmeans-n-clusters <number>                     Target at <number> clusters for K means. [default: 20]

  --run-hdbscan                                    Run hdbscan clustering algorithm on diffusion components.
  --hdbscan-min-cluster-size <number>              Minimum cluster size for hdbscan. [default: 50]
  --hdbscan-min-samples <number>                   Minimum number of samples for hdbscan. [default: 50]

  --run-approximated-louvain                       Run approximated louvain clustering algorithm.
  --approx-louvain-ninit <number>                  Number of Kmeans tries. [default: 20]
  --approx-louvain-nclusters <number>              Number of clusters for Kmeans initialization. [default: 30]
  --approx-louvain-resolution <resolution>.        Resolution parameter for louvain. [default: 1.3]

  --run-tsne                                       Run multi-core t-SNE for visualization.
  --tsne-perplexity <perplexity>                   t-SNE's perplexity parameter. [default: 30]
  --run-fitsne                                     Run FIt-SNE for visualization.

  --run-umap                                       Run umap for visualization.
  --umap-on-diffmap                                Run umap on diffusion components.
  --umap-K <K>                                     K neighbors for umap. [default: 15]
  --umap-min-dist <number>                         Umap parameter. [default: 0.1]
  --umap-spread <spread>                           Umap parameter. [default: 1.0]

  --run-fle                                        Run force-directed layout embedding.
  --fle-K <K>                                      K neighbors for building graph for FLE. [default: 50]
  --fle-n-steps <nstep>                            Number of iterations for FLE. [default: 10000]
  --fle-affinity <affinity>                        Affinity matrix to be used. Could be 'W_diffmap', or 'W_diffmap_norm'. [default: W_diffmap]

  -h, --help                                       Print out help information.

Outputs:
  output_name.h5ad              Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, t-SNE and diffusion map coordinates are stored in the 'obsm' field.
  output_name.seurat.h5ad       Optional output. Only exists if '--output-seurat-compatible' is set. This is the Seurat-readable h5ad file.
  output_name.loom              Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.
  
Examples:
  scCloud cluster -p 20 --correct-batch-effect --run-louvain --run-tsne manton_bm_10x.h5 manton_bm
    """

    def execute(self):
        kwargs = {
            'n_jobs' : int(self.args['--threads']),
            'genome' : self.args['--genome'],

            'processed' : self.args['--processed'],
            'subcluster' : False,

            'cite_seq' : self.args['--cite-seq'],

            'filt_xlsx' : self.args['--output-filtration-results'],
            'output_seurat_compatible' : self.args['--output-seurat-compatible'],
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
            'diffmap_full_speed' : self.args['--diffmap-full-speed'],

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
