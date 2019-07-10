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

  --select-hvg-flavor <flavor>                     Highly variable gene selection method. <flavor> can be 'scCloud' or 'Seurat'. [default: scCloud]
  --select-hvg-ngenes <ngenes>                     Select top <ngenes> highly variable genes. If <flavor> is 'Seurat' and <ngenes> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]
  --no-select-hvg                                  Do not select highly variable genes.
  --plot-hvg                                       Plot highly variable gene selection.

  --random-state <seed>                            Random number generator seed. [default: 0]
  --temp-folder <temp_folder>                      Joblib temporary folder for memmapping numpy arrays.

  --nPC <number>                                   Number of PCs. [default: 50]

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
  --approx-leiden-class-label <label>              Approximated leiden label name in AnnData. [default: approx_louvain_labels]

  --run-tsne                                       Run multi-core t-SNE for visualization.
  --run-fitsne                                     Run FIt-SNE for visualization.
  --tsne-perplexity <perplexity>                   t-SNE's perplexity parameter, used by both tSNE, FItSNE net-tSNE and net-FItSNE. [default: 30]

  --run-umap                                       Run umap for visualization.
  --umap-K <K>                                     K neighbors for umap. [default: 15]
  --umap-min-dist <number>                         Umap parameter. [default: 0.1]
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
            'select_singlets' : False,
            
            'subset_selections' : self.args['--subset-selection'],

            'n_jobs' : int(self.args['--threads']),
            'genome' : None,
            'batch_correction' : self.args['--correct-batch-effect'],
            'output_loom' : self.args['--output-loom'],

            'select_hvg' : not self.args['--no-select-hvg'],
            'hvg_flavor' : self.args['--select-hvg-flavor'],
            'hvg_ngenes' : int(self.args['--select-hvg-ngenes']) if self.args['--select-hvg-ngenes'] != 'None' else None,
            'plot_hvg' : self.args['<output_name>'] if self.args['--plot-hvg'] else None,

            'random_state' : int(self.args['--random-state']),
            'temp_folder' : self.args['--temp-folder'],

            'nPC' : int(self.args['--nPC']),

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
