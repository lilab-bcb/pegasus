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

  --select-singlets                                Only select DemuxEM-predicted singlets for analysis.  
  --cite-seq                                       Data are CITE-Seq data. scCloud will perform analyses on RNA count matrix first. Then it will attach the ADT matrix to the RNA matrix with all antibody names changing to 'AD-' + antibody_name. Lastly, it will embed the antibody expression using FIt-SNE (the basis used for plotting is 'citeseq_fitsne').
  --cite-seq-capping <percentile>                  For CITE-Seq surface protein expression, make all cells with expression > <percentile> to the value at <percentile> to smooth outlier. Set <percentile> to 100.0 to turn this option off. [default: 99.99]

  --output-filtration-results                      Output filtration results as a spreadsheet.
  --plot-filtration-results                        Plot filtration results as PDF files.
  --plot-filtration-figsize <figsize>              Figure size for filtration plots. <figsize> is a comma-separated list of two numbers, the width and height of the figure (e.g. 6,4).
  --output-seurat-compatible                       Output seurat-compatible h5ad file. Caution: File size might be large, do not turn this option on for large data sets.
  --output-loom                                    Output loom-formatted file.
  
  --min-genes <number>                             Only keep cells with at least <number> of genes. [default: 500]
  --max-genes <number>                             Only keep cells with less than <number> of genes. [default: 6000]
  --min-umis <number>                              Only keep cells with at least <number> of UMIs. [default: 100]
  --max-umis <number>                              Only keep cells with less than <number> of UMIs. [default: 600000]
  --mito-prefix <prefix>                           Prefix for mitochondrial genes. If multiple prefixes are provided, separate them by comma (e.g. "MT-,mt-"). [default: MT-]
  --percent-mito <ratio>                           Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]
  --gene-percent-cells <ratio>                     Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]
  --min-genes-on-raw <number>                      If input are raw 10x matrix, which include all barcodes, perform a pre-filtration step to keep the data size small. In the pre-filtration step, only keep cells with at least <number> of genes. [default: 100]

  --counts-per-cell-after <number>                 Total counts per cell after normalization. [default: 1e5]

  --correct-batch-effect                           Correct for batch effects.
  --batch-group-by <expression>                    Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+...+attrn', or 'attr=value11,...,value1n_1;value21,...,value2n_2;...;valuem1,...,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',...,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,...,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.

  --select-hvg-flavor <flavor>                     Highly variable gene selection method. <flavor> can be 'scCloud' or 'Seurat'. [default: scCloud]
  --select-hvg-ngenes <ngenes>                     Select top <ngenes> highly variable genes. If <flavor> is 'Seurat' and <ngenes> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]
  --benchmark-time                                 This option is used for benchmarking time, will calculate mean and variance even if they are calculated in batch correction.
  --no-select-hvg                                  Do not select highly variable genes.
  --plot-hvg                                       Plot highly variable gene selection.

  --random-state <seed>                            Random number generator seed. [default: 0]
  --temp-folder <temp_folder>                      Joblib temporary folder for memmapping numpy arrays.

  --nPC <number>                                   Number of principal components. [default: 50]

  --kBET                                           Calculate kBET.
  --kBET-batch <batch>                             kBET batch keyword.
  --kBET-alpha <alpha>                             kBET rejection alpha. [default: 0.05]
  --kBET-K <K>                                     kBET K. [default: 25]

  --nDC <number>                                   Number of diffusion components. [default: 50]
  --diffmap-alpha <alpha>                          Power parameter for diffusion-based pseudotime. [default: 0.5]
  --diffmap-K <K>                                  Number of neighbors used for constructing affinity matrix. [default: 100]
  --diffmap-full-speed                             For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.
  --diffmap-solver <solver>                        Solver for eigen decomposition, either 'randomized' or 'eigsh'. [default: randomized]

  --calculate-pseudotime <roots>                   Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  --run-louvain                                    Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --louvain-affinity <affinity>                    Affinity matrix to be used. Could be 'W' or 'W_diffmap'. [default: W]
  --louvain-class-label <label>                    Louvain cluster label name in AnnData. [default: louvain_labels]

  --run-leiden                                     Run leiden clustering algorithm.
  --leiden-resolution <resolution>                 Resolution parameter for the leiden clustering algorithm. [default: 1.3]
  --leiden-affinity <affinity>                     Affinity matrix to be used. Could be 'W' or 'W_diffmap'. [default: W]
  --leiden-niter <niter>                           Number of iterations of running the Leiden algorithm. If <niter> is negative, run Leiden iteratively until no improvement. [default: -1]
  --leiden-class-label <label>                     Leiden cluster label name in AnnData. [default: leiden_labels]

  --run-approximated-louvain                       Run approximated louvain clustering algorithm.
  --approx-louvain-basis <basis>                   Basis used for KMeans clustering. Can be 'pca', 'rpca', or 'diffmap'. [default: diffmap]
  --approx-louvain-nclusters <number>              Number of clusters for Kmeans initialization. [default: 30]
  --approx-louvain-ninit <number>                  Number of Kmeans tries. [default: 20]
  --approx-louvain-resolution <resolution>         Resolution parameter for louvain. [default: 1.3]
  --approx-louvain-affinity <affinity>             Affinity matrix to be used. Could be 'W' or 'W_diffmap'. [default: W]
  --approx-louvain-class-label <label>             Approximated louvain label name in AnnData. [default: approx_louvain_labels]

  --run-approximated-leiden                        Run approximated leiden clustering algorithm.
  --approx-leiden-basis <basis>                    Basis used for KMeans clustering. Can be 'pca', 'rpca', or 'diffmap'. [default: diffmap]
  --approx-leiden-nclusters <number>               Number of clusters for Kmeans initialization. [default: 30]
  --approx-leiden-ninit <number>                   Number of Kmeans tries. [default: 20]
  --approx-leiden-resolution <resolution>          Resolution parameter for leiden. [default: 1.3]
  --approx-leiden-affinity <affinity>              Affinity matrix to be used. Could be 'W' or 'W_diffmap'. [default: W]
  --approx-leiden-class-label <label>              Approximated leiden label name in AnnData. [default: approx_leiden_labels]

  --run-tsne                                       Run multi-core t-SNE for visualization.
  --run-fitsne                                     Run FIt-SNE for visualization.
  --tsne-perplexity <perplexity>                   t-SNE's perplexity parameter, used by both tSNE, FItSNE net-tSNE and net-FItSNE. [default: 30]

  --run-umap                                       Run umap for visualization.
  --umap-K <K>                                     K neighbors for umap. [default: 15]
  --umap-min-dist <number>                         Umap parameter. [default: 0.5]
  --umap-spread <spread>                           Umap parameter. [default: 1.0]

  --run-fle                                        Run force-directed layout embedding.
  --fle-K <K>                                      K neighbors for building graph for FLE. [default: 50]
  --fle-target-change-per-node <change>            Target change per node to stop forceAtlas2. [default: 2.0]
  --fle-target-steps <steps>                       Maximum number of iterations before stopping the forceAtlas2 algoritm. [default: 5000]
  --fle-3D                                         Calculate 3D force-directed layout. 

  --net-down-sample-fraction <frac>                Down sampling fraction for net-related visualization. [default: 0.1]
  --net-down-sample-K <K>                          Use <K> neighbors to estimate local density for each data point for down sampling. [default: 25]
  --net-down-sample-alpha <alpha>                  Weighted down sample, proportional to radius^alpha. [default: 1.0]

  --net-regressor-L2-penalty <value>               L2 penalty parameter for the deep net regressor. [default: 0.1]
  --net-ds-full-speed                              For net-UMAP and net-FLE, use full speed for the down-sampled data.

  --run-net-tsne                                   Run net tSNE for visualization.
  --net-tsne-polish-learning-frac <frac>           After running the deep regressor to predict new coordinates, use <frac> * nsample as the learning rate to use to polish the coordinates. [default: 0.33]
  --net-tsne-polish-niter <niter>                  Number of iterations for polishing tSNE run. [default: 150]
  --net-tsne-out-basis <basis>                     Output basis for net-tSNE. [default: net_tsne]

  --run-net-fitsne                                 Run net FIt-SNE for visualization.
  --net-fitsne-polish-learning-frac <frac>         After running the deep regressor to predict new coordinates, use <frac> * nsample as the learning rate to use to polish the coordinates. [default: 0.5]
  --net-fitsne-polish-niter <niter>                Number of iterations for polishing FItSNE run. [default: 150]
  --net-fitsne-out-basis <basis>                   Output basis for net-FItSNE. [default: net_fitsne]

  --run-net-umap                                   Run net umap for visualization.
  --net-umap-polish-learning-rate <rate>           After running the deep regressor to predict new coordinate, what is the learning rate to use to polish the coordinates for UMAP. [default: 1.0]
  --net-umap-polish-nepochs <nepochs>              Number of iterations for polishing UMAP run. [default: 40]
  --net-umap-out-basis <basis>                     Output basis for net-UMAP. [default: net_umap]

  --run-net-fle                                    Run net FLE.
  --net-fle-ds-full-speed                          If run full-speed kNN on down-sampled data points.
  --net-fle-polish-target-steps <steps>            After running the deep regressor to predict new coordinate, what is the number of force atlas 2 iterations. [default: 1500]
  --net-fle-out-basis <basis>                      Output basis for net-FLE. [default: net_fle]

  -h, --help                                       Print out help information.

Outputs:
  output_name.h5ad                 Output file in h5ad format. To load this file in python, use ``import scCloud; data = scCloud.tools.read_input('output_name.h5ad', mode = 'a')``. The log-normalized expression matrix is stored in ``data.X`` as a CSR-format sparse matrix. The ``obs`` field contains cell related attributes, including clustering results. For example, ``data.obs_names`` records cell barcodes; ``data.obs['Channel']`` records the channel each cell comes from; ``data.obs['n_genes']``, ``data.obs['n_counts']``, and ``data.obs['percent_mito']`` record the number of expressed genes, total UMI count, and mitochondrial rate for each cell respectively; ``data.obs['louvain_labels']`` and ``data.obs['approx_louvain_labels']`` record each cell's cluster labels using different clustring algorithms; ``data.obs['pseudo_time']`` records the inferred pseudotime for each cell. The ``var`` field contains gene related attributes. For example, ``data.var_names`` records gene symbols, ``data.var['gene_ids']`` records Ensembl gene IDs, and ``data.var['selected']`` records selected variable genes. The ``obsm`` field records embedding coordiates. For example, ``data.obsm['X_pca']`` records PCA coordinates, ``data.obsm['X_tsne']`` records tSNE coordinates, ``data.obsm['X_umap']`` records UMAP coordinates, ``data.obsm['X_diffmap']`` records diffusion map coordinates, ``data.obsm['X_diffmap_pca']`` records the first 3 PCs by projecting the diffusion components using PCA, and ``data.obsm['X_fle']`` records the force-directed layout coordinates from the diffusion components. The ``uns`` field stores other related information, such as reference genome (``data.uns['genome']``). If '--make-output-seurat-compatible' is on, this file can be loaded into R and converted into a Seurat object.
  output_name.seurat.h5ad          Optional output. Only exists if '--output-seurat-compatible' is set. 'output_name.h5ad' in seurat-compatible manner. This file can be loaded into R and converted into a Seurat object.
  output_name.filt.xlsx            Optional output. Only exists if '--output-filtration-results' is set. This file has two sheets --- Cell filtration stats and Gene filtration stats. The first sheet records cell filtering results and it has 10 columns: Channel, channel name; kept, number of cells kept; median_n_genes, median number of expressed genes in kept cells; median_n_umis, median number of UMIs in kept cells; median_percent_mito, median mitochondrial rate as UMIs between mitochondrial genes and all genes in kept cells; filt, number of cells filtered out; total, total number of cells before filtration, if the input contain all barcodes, this number is the cells left after '--min-genes-on-raw' filtration; median_n_genes_before, median expressed genes per cell before filtration; median_n_umis_before, median UMIs per cell before filtration; median_percent_mito_before, median mitochondrial rate per cell before filtration. The channels are sorted in ascending order with respect to the number of kept cells per channel. The second sheet records genes that failed to pass the filtering. This sheet has 3 columns: gene, gene name; n_cells, number of cells this gene is expressed; percent_cells, the fraction of cells this gene is expressed. Genes are ranked in ascending order according to number of cells the gene is expressed. Note that only genes not expressed in any cell are removed from the data. Other filtered genes are marked as non-robust and not used for TPM-like normalization.
  output_name.filt.gene.pdf        Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting gene count distributions before and after filtration per channel.
  output_name.filt.UMI.pdf         Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting UMI count distributions before and after filtration per channel.
  output_name.filt.mito.pdf        Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting mitochondrial rate distributions before and after filtration per channel.
  output_name.hvg.pdf              Optional output. Only exists if '--plot-hvg' is set. This file contains a scatter plot describing the highly variable gene selection procedure.
  output_name.loom                 Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.
  
Examples:
  scCloud cluster -p 20 --correct-batch-effect --run-louvain --run-tsne manton_bm_10x.h5 manton_bm
    """

    def execute(self):
        kwargs = {
            'n_jobs' : int(self.args['--threads']),
            'genome' : self.args['--genome'],

            'processed' : self.args['--processed'],
            'subcluster' : False,

            'select_singlets' : self.args['--select-singlets'],
            'cite_seq' : self.args['--cite-seq'],
            'cite_seq_capping' : float(self.args['--cite-seq-capping']),

            'output_filt' : self.args['<output_name>'] if self.args['--output-filtration-results'] else None,
            'plot_filt' : self.args['<output_name>'] if self.args['--plot-filtration-results'] else None,
            'plot_filt_figsize' : self.args['--plot-filtration-figsize'],

            'seurat_compatible' : self.args['--output-seurat-compatible'],
            'output_loom' : self.args['--output-loom'],

            'min_genes' : int(self.args['--min-genes']),
            'max_genes' : int(self.args['--max-genes']),
            'min_umis' : int(self.args['--min-umis']),
            'max_umis' : int(self.args['--max-umis']),
            'mito_prefix' : self.args['--mito-prefix'],
            'percent_mito' : float(self.args['--percent-mito']),
            'percent_cells' : float(self.args['--gene-percent-cells']),
            'min_genes_on_raw' : int(self.args['--min-genes-on-raw']),

            'norm_count' : float(self.args['--counts-per-cell-after']),

            'batch_correction' : self.args['--correct-batch-effect'],
            'group_attribute' : self.args['--batch-group-by'],

            'select_hvg' : not self.args['--no-select-hvg'],
            'hvg_flavor' : self.args['--select-hvg-flavor'],
            'hvg_ngenes' : int(self.args['--select-hvg-ngenes']) if self.args['--select-hvg-ngenes'] != 'None' else None,
            'benchmark_time' : self.args['--benchmark-time'],
            'plot_hvg' : self.args['<output_name>'] if self.args['--plot-hvg'] else None,

            'random_state' : int(self.args['--random-state']),
            'temp_folder' : self.args['--temp-folder'],
            
            'nPC' : int(self.args['--nPC']),

            'kBET' : self.args['--kBET'],
            'kBET_batch' : self.args['--kBET-batch'],
            'kBET_alpha' : float(self.args['--kBET-alpha']),
            'kBET_K': int(self.args['--kBET-K']),

            'nDC' : int(self.args['--nDC']),
            'diffmap_alpha' : float(self.args['--diffmap-alpha']),
            'diffmap_K' : int(self.args['--diffmap-K']),
            'diffmap_full_speed' : self.args['--diffmap-full-speed'],
            'diffmap_solver' : self.args['--diffmap-solver'],

            'run_louvain' : self.args['--run-louvain'],
            'louvain_resolution' : float(self.args['--louvain-resolution']),
            'louvain_affinity' : self.args['--louvain-affinity'],
            'louvain_class_label' : self.args['--louvain-class-label'],

            'run_leiden' : self.args['--run-leiden'],
            'leiden_resolution' : float(self.args['--leiden-resolution']),
            'leiden_affinity' : self.args['--leiden-affinity'],
            'leiden_niter' : int(self.args['--leiden-niter']),
            'leiden_class_label' : self.args['--leiden-class-label'],

            'run_approx_louvain' : self.args['--run-approximated-louvain'],
            'approx_louvain_basis' : self.args['--approx-louvain-basis'],
            'approx_louvain_nclusters' : int(self.args['--approx-louvain-nclusters']),
            'approx_louvain_ninit' : int(self.args['--approx-louvain-ninit']),
            'approx_louvain_resolution' : float(self.args['--approx-louvain-resolution']),
            'approx_louvain_affinity' : self.args['--approx-louvain-affinity'],
            'approx_louvain_class_label' : self.args['--approx-louvain-class-label'],

            'run_approx_leiden' : self.args['--run-approximated-leiden'],
            'approx_leiden_basis' : self.args['--approx-leiden-basis'],
            'approx_leiden_nclusters' : int(self.args['--approx-leiden-nclusters']),
            'approx_leiden_ninit' : int(self.args['--approx-leiden-ninit']),
            'approx_leiden_resolution' : float(self.args['--approx-leiden-resolution']),
            'approx_leiden_affinity' : self.args['--approx-leiden-affinity'],
            'approx_leiden_class_label' : self.args['--approx-leiden-class-label'],

            'run_tsne' : self.args['--run-tsne'],
            'run_fitsne' : self.args['--run-fitsne'],
            'tsne_perplexity' : float(self.args['--tsne-perplexity']),

            'run_umap' : self.args['--run-umap'],
            'umap_K' : int(self.args['--umap-K']),
            'umap_min_dist' : float(self.args['--umap-min-dist']),
            'umap_spread' : float(self.args['--umap-spread']),

            'run_fle' : self.args['--run-fle'],
            'fle_K' : int(self.args['--fle-K']),
            'fle_target_change_per_node' : float(self.args['--fle-target-change-per-node']),
            'fle_target_steps' : int(self.args['--fle-target-steps']),
            'fle_3D' : self.args['--fle-3D'],

            'net_ds_frac' : float(self.args['--net-down-sample-fraction']),
            'net_ds_K' : int(self.args['--net-down-sample-K']),
            'net_ds_alpha' : float(self.args['--net-down-sample-alpha']),

            'net_l2' : float(self.args['--net-regressor-L2-penalty']),
            'net_ds_full_speed' : self.args['--net-ds-full-speed'],

            'run_net_tsne' : self.args['--run-net-tsne'],
            'net_tsne_polish_learing_frac' : float(self.args['--net-tsne-polish-learning-frac']),
            'net_tsne_polish_niter' : int(self.args['--net-tsne-polish-niter']),
            'net_tsne_basis' : self.args['--net-tsne-out-basis'],

            'run_net_fitsne' : self.args['--run-net-fitsne'],
            'net_fitsne_polish_learing_frac' : float(self.args['--net-fitsne-polish-learning-frac']),
            'net_fitsne_polish_niter' : int(self.args['--net-fitsne-polish-niter']),
            'net_fitsne_basis' : self.args['--net-fitsne-out-basis'],

            'run_net_umap' : self.args['--run-net-umap'],
            'net_umap_polish_learing_rate' : float(self.args['--net-umap-polish-learning-rate']),
            'net_umap_polish_nepochs' : int(self.args['--net-umap-polish-nepochs']),
            'net_umap_basis' : self.args['--net-umap-out-basis'],

            'run_net_fle' : self.args['--run-net-fle'],
            'net_fle_ds_full_speed' : self.args['--net-fle-ds-full-speed'],
            'net_fle_polish_target_steps' : int(self.args['--net-fle-polish-target-steps']),
            'net_fle_basis' : self.args['--net-fle-out-basis'],

            'pseudotime' : self.split_string(self.args['--calculate-pseudotime'])
        }

        run_pipeline(self.args['<input_file>'], self.args['<output_name>'], **kwargs)
