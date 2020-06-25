from .Base import Base
from pegasus.pipeline import run_pipeline


class SubClustering(Base):
    """
Run pegasus to obtain subclusters.

Usage:
  pegasus subcluster [options] --subset-selection <subset-selection>... <input_file> <output_name>
  pegasus subcluster -h

Arguments:
  input_file             Single cell data with clustering done in h5ad format.
  output_name            Output file name. All outputs will use it as the prefix.

Options:
  --subset-selection <subset-selection>...         Specify which cells will be included in the subcluster analysis. Each <subset_selection> string takes the format of 'attr:value,...,value', which means select cells with attr in the values. If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings.

  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --batch-group-by                                 Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+...+attrn', or 'attr=value11,...,value1n_1;value21,...,value2n_2;...;valuem1,...,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',...,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,...,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.

  --output-loom                                    Output loom-formatted file.

  --select-hvf-flavor <flavor>                     Highly variable feature selection method. <flavor> can be 'pegasus' or 'Seurat'. [default: pegasus]
  --select-hvf-ngenes <nfeatures>                  Select top <nfeatures> highly variable features. If <flavor> is 'Seurat' and <nfeatures> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]
  --no-select-hvf                                  Do not select highly variable features.
  --plot-hvf                                       Plot highly variable feature selection.

  --correct-batch-effect                           Correct for batch effects for subclustering task.
  --correction-method <method>                     Batch correction method, can be either 'L/S' for location/scale adjustment algorithm (Li and Wong. The analysis of Gene Expression Data 2003) or 'harmony' for Harmony (Korsunsky et al. Nature Methods 2019) or 'scanorama' for Scanorama (Hie et al. Nature Biotechnology 2019). For L/S method, we directly use correction factors from the parent study. [default: harmony]
  --harmony-nclusters <nclusters>                  Number of clusters used for Harmony batch correction.

  --random-state <seed>                            Random number generator seed. [default: 0]
  --temp-folder <temp_folder>                      Joblib temporary folder for memmapping numpy arrays.

  --pca-n <number>                                 Number of principal components. [default: 50]
  --pca-robust                                     Use 'arpack' instead of 'randomized' as svd_solver for large sparse matrices. It will take longer time to compute PCs, but the results are more numerically stable. 

  --knn-K <number>                                 Number of nearest neighbors for building kNN graph. [default: 100]
  --knn-full-speed                                 For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

  --kBET                                           Calculate kBET.
  --kBET-batch <batch>                             kBET batch keyword.
  --kBET-alpha <alpha>                             kBET rejection alpha. [default: 0.05]
  --kBET-K <K>                                     kBET K. [default: 25]

  --diffmap                                        Calculate diffusion maps.
  --diffmap-ndc <number>                           Number of diffusion components. [default: 100]
  --diffmap-solver <solver>                        Solver for eigen decomposition, either 'eigsh' or 'randomized'. [default: eigsh]
  --diffmap-maxt <max_t>                           Maximum time stamp to search for the knee point. [default: 5000]
  --diffmap-to-3d                                  If map diffusion map into 3D space using PCA.
  --calculate-pseudotime <roots>                   Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  --louvain                                        Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --louvain-class-label <label>                    Louvain cluster label name in AnnData. [default: louvain_labels]

  --leiden                                         Run leiden clustering algorithm.
  --leiden-resolution <resolution>                 Resolution parameter for the leiden clustering algorithm. [default: 1.3]
  --leiden-niter <niter>                           Number of iterations of running the Leiden algorithm. If <niter> is negative, run Leiden iteratively until no improvement. [default: -1]
  --leiden-class-label <label>                     Leiden cluster label name in AnnData. [default: leiden_labels]

  --spectral-louvain                               Run spectral-louvain clustering algorithm.
  --spectral-louvain-basis <basis>                 Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. If 'diffmap' is not calculated, use 'pca' instead. [default: diffmap]
  --spectral-louvain-nclusters <number>            Number of first level clusters for Kmeans. [default: 30]
  --spectral-louvain-nclusters2 <number>           Number of second level clusters for Kmeans. [default: 50]
  --spectral-louvain-ninit <number>                Number of Kmeans tries for first level clustering. Default is the same as scikit-learn Kmeans function. [default: 10]
  --spectral-louvain-resolution <resolution>       Resolution parameter for louvain. [default: 1.3]
  --spectral-louvain-class-label <label>           Spectral-louvain label name in AnnData. [default: spectral_louvain_labels]

  --spectral-leiden                                Run spectral-leiden clustering algorithm.
  --spectral-leiden-basis <basis>                  Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. If 'diffmap' is not calculated, use 'pca' instead. [default: diffmap]
  --spectral-leiden-nclusters <number>             Number of first level clusters for Kmeans. [default: 30]
  --spectral-leiden-nclusters2 <number>            Number of second level clusters for Kmeans. [default: 50]
  --spectral-leiden-ninit <number>                 Number of Kmeans tries for first level clustering. Default is the same as scikit-learn Kmeans function. [default: 10]
  --spectral-leiden-resolution <resolution>        Resolution parameter for leiden. [default: 1.3]
  --spectral-leiden-class-label <label>            Spectral-leiden label name in AnnData. [default: spectral_leiden_labels]

  --tsne                                           Run multi-core t-SNE for visualization.
  --fitsne                                         Run FIt-SNE for visualization.
  --tsne-perplexity <perplexity>                   t-SNE's perplexity parameter, used by both tSNE, FItSNE net-tSNE and net-FItSNE. [default: 30]

  --umap                                           Run umap for visualization.
  --umap-K <K>                                     K neighbors for umap. [default: 15]
  --umap-min-dist <number>                         Umap parameter. [default: 0.5]
  --umap-spread <spread>                           Umap parameter. [default: 1.0]

  --fle                                            Run force-directed layout embedding.
  --fle-K <K>                                      K neighbors for building graph for FLE. [default: 50]
  --fle-target-change-per-node <change>            Target change per node to stop forceAtlas2. [default: 2.0]
  --fle-target-steps <steps>                       Maximum number of iterations before stopping the forceAtlas2 algoritm. [default: 5000]
  --fle-memory <memory>                            Memory size in GB for the Java FA2 component. [default: 8]

  --net-down-sample-fraction <frac>                Down sampling fraction for net-related visualization. [default: 0.1]
  --net-down-sample-K <K>                          Use <K> neighbors to estimate local density for each data point for down sampling. [default: 25]
  --net-down-sample-alpha <alpha>                  Weighted down sample, proportional to radius^alpha. [default: 1.0]

  --net-regressor-L2-penalty <value>               L2 penalty parameter for the deep net regressor. [default: 0.1]

  --net-tsne                                       Run net tSNE for visualization.
  --net-tsne-polish-learning-frac <frac>           After running the deep regressor to predict new coordinates, use <frac> * nsample as the learning rate to use to polish the coordinates. [default: 0.33]
  --net-tsne-polish-niter <niter>                  Number of iterations for polishing tSNE run. [default: 150]
  --net-tsne-out-basis <basis>                     Output basis for net-tSNE. [default: net_tsne]

  --net-umap                                       Run net umap for visualization.
  --net-umap-polish-learning-rate <rate>           After running the deep regressor to predict new coordinate, what is the learning rate to use to polish the coordinates for UMAP. [default: 1.0]
  --net-umap-polish-nepochs <nepochs>              Number of iterations for polishing UMAP run. [default: 40]
  --net-umap-out-basis <basis>                     Output basis for net-UMAP. [default: net_umap]

  --net-fle                                        Run net FLE.
  --net-fle-polish-target-steps <steps>            After running the deep regressor to predict new coordinate, what is the number of force atlas 2 iterations. [default: 1500]
  --net-fle-out-basis <basis>                      Output basis for net-FLE. [default: net_fle]

  -h, --help                                       Print out help information.

Outputs:
  output_name.h5ad              Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, t-SNE and diffusion map coordinates are stored in the 'obsm' field.
  output_name.loom              Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.

Examples:
  pegasus subcluster -p 20 --correct-batch-effect --subset-selection louvain_labels:3,6 --subset-selection Condition:CB_nonmix --tsne --louvain manton_bm.h5ad manton_bm_subset
    """

    def execute(self):
        kwargs = {
            "processed": True,
            "subcluster": True,
            "cite_seq": False,
            "select_singlets": False,
            "subset_selections": self.args["--subset-selection"],
            "n_jobs": int(self.args["--threads"]),
            "genome": None,
            "channel_attr": None,
            "black_list": None,
            "batch_correction": self.args["--correct-batch-effect"],
            "group_attribute": self.args["--batch-group-by"],
            "output_loom": self.args["--output-loom"],
            "select_hvf": not self.args["--no-select-hvf"],
            "hvf_flavor": self.args["--select-hvf-flavor"],
            "hvf_ngenes": int(self.args["--select-hvf-ngenes"])
            if self.args["--select-hvf-ngenes"] != "None"
            else None,
            "plot_hvf": self.args["<output_name>"] if self.args["--plot-hvf"] else None,
            "batch_correction": self.args["--correct-batch-effect"],
            "correction_method": self.args["--correction-method"],
            "harmony_nclusters": self.convert_to_int(self.args["--harmony-nclusters"]),
            "random_state": int(self.args["--random-state"]),
            "temp_folder": self.args["--temp-folder"],
            "pca_n": int(self.args["--pca-n"]),
            "pca_robust": self.args["--pca-robust"],
            "K": int(self.args["--knn-K"]),
            "full_speed": self.args["--knn-full-speed"],
            "kBET": self.args["--kBET"],
            "kBET_batch": self.args["--kBET-batch"],
            "kBET_alpha": float(self.args["--kBET-alpha"]),
            "kBET_K": int(self.args["--kBET-K"]),
            "diffmap": self.args["--diffmap"],
            "diffmap_ndc": int(self.args["--diffmap-ndc"]),
            "diffmap_maxt": int(self.args["--diffmap-maxt"]),
            "diffmap_solver": self.args["--diffmap-solver"],
            "diffmap_to_3d": self.args["--diffmap-to-3d"],
            "pseudotime": self.split_string(self.args["--calculate-pseudotime"]),
            "louvain": self.args["--louvain"],
            "louvain_resolution": float(self.args["--louvain-resolution"]),
            "louvain_class_label": self.args["--louvain-class-label"],
            "leiden": self.args["--leiden"],
            "leiden_resolution": float(self.args["--leiden-resolution"]),
            "leiden_niter": int(self.args["--leiden-niter"]),
            "leiden_class_label": self.args["--leiden-class-label"],
            "spectral_louvain": self.args["--spectral-louvain"],
            "spectral_louvain_basis": self.args["--spectral-louvain-basis"],
            "spectral_louvain_nclusters": int(
                self.args["--spectral-louvain-nclusters"]
            ),
            "spectral_louvain_nclusters2": int(
                self.args["--spectral-louvain-nclusters2"]
            ),
            "spectral_louvain_ninit": int(self.args["--spectral-louvain-ninit"]),
            "spectral_louvain_resolution": float(
                self.args["--spectral-louvain-resolution"]
            ),
            "spectral_louvain_class_label": self.args["--spectral-louvain-class-label"],
            "spectral_leiden": self.args["--spectral-leiden"],
            "spectral_leiden_basis": self.args["--spectral-leiden-basis"],
            "spectral_leiden_nclusters": int(self.args["--spectral-leiden-nclusters"]),
            "spectral_leiden_nclusters2": int(self.args["--spectral-leiden-nclusters2"]),
            "spectral_leiden_ninit": int(self.args["--spectral-leiden-ninit"]),
            "spectral_leiden_resolution": float(
                self.args["--spectral-leiden-resolution"]
            ),
            "spectral_leiden_class_label": self.args["--spectral-leiden-class-label"],
            "tsne": self.args["--tsne"],
            "fitsne": self.args["--fitsne"],
            "tsne_perplexity": float(self.args["--tsne-perplexity"]),
            "umap": self.args["--umap"],
            "umap_K": int(self.args["--umap-K"]),
            "umap_min_dist": float(self.args["--umap-min-dist"]),
            "umap_spread": float(self.args["--umap-spread"]),
            "fle": self.args["--fle"],
            "fle_K": int(self.args["--fle-K"]),
            "fle_target_change_per_node": float(
                self.args["--fle-target-change-per-node"]
            ),
            "fle_target_steps": int(self.args["--fle-target-steps"]),
            "fle_memory": int(self.args["--fle-memory"]),
            "net_ds_frac": float(self.args["--net-down-sample-fraction"]),
            "net_ds_K": int(self.args["--net-down-sample-K"]),
            "net_ds_alpha": float(self.args["--net-down-sample-alpha"]),
            "net_l2": float(self.args["--net-regressor-L2-penalty"]),
            "net_tsne": self.args["--net-tsne"],
            "net_tsne_polish_learing_frac": float(
                self.args["--net-tsne-polish-learning-frac"]
            ),
            "net_tsne_polish_niter": int(self.args["--net-tsne-polish-niter"]),
            "net_tsne_basis": self.args["--net-tsne-out-basis"],
            "net_umap": self.args["--net-umap"],
            "net_umap_polish_learing_rate": float(
                self.args["--net-umap-polish-learning-rate"]
            ),
            "net_umap_polish_nepochs": int(self.args["--net-umap-polish-nepochs"]),
            "net_umap_basis": self.args["--net-umap-out-basis"],
            "net_fle": self.args["--net-fle"],
            "net_fle_polish_target_steps": int(
                self.args["--net-fle-polish-target-steps"]
            ),
            "net_fle_basis": self.args["--net-fle-out-basis"],
        }

        import logging
        logger = logging.getLogger("pegasus")
        fh = logging.FileHandler("{}.log".format(self.args["<output_name>"]))
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        run_pipeline(self.args["<input_file>"], self.args["<output_name>"], **kwargs)
