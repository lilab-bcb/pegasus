from .Base import Base
from pegasus.pipeline import run_pipeline


class Clustering(Base):
    """
Run pegasus.pipeline to obtain top-level clusters.

Usage:
  pegasus cluster [options] <input_file> <output_name>
  pegasus cluster -h

Arguments:
  input_file       Input file in either 'zarr', 'h5ad', 'loom', '10x', 'mtx', 'csv', 'tsv' or 'fcs' format. If first-pass analysis has been performed, but you want to run some additional analysis, you could also pass a zarr-formatted file.
  output_name      Output file name. All outputs will use it as the prefix.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --processed                                      Input file is processed. Assume quality control, data normalization and log transformation, highly variable gene selection, batch correction/PCA and kNN graph building is done.

  --black-list <black_list>                        Cell barcode attributes in black list will be popped out. Format is "attr1,attr2,...,attrn".

  --select-singlets                                Only select DemuxEM-predicted singlets for analysis.
  --remap-singlets <remap_string>                  Remap singlet names using <remap_string>, where <remap_string> takes the format "new_name_i:old_name_1,old_name_2;new_name_ii:old_name_3;...". For example, if we hashed 5 libraries from 3 samples sample1_lib1, sample1_lib2, sample2_lib1, sample2_lib2 and sample3, we can remap them to 3 samples using this string: "sample1:sample1_lib1,sample1_lib2;sample2:sample2_lib1,sample2_lib2". In this way, the new singlet names will be in metadata field with key 'assignment', while the old names will be kept in metadata field with key 'assignment.orig'.
  --subset-singlets <subset_string>                If select singlets, only select singlets in the <subset_string>, which takes the format "name1,name2,...". Note that if --remap-singlets is specified, subsetting happens after remapping. For example, we can only select singlets from sampe 1 and 3 using "sample1,sample3".

  --genome <genome_name>                           If sample count matrix is in either DGE, mtx, csv, tsv or loom format, use <genome_name> as the genome reference name.
  --focus <keys>                                   Focus analysis on Unimodal data with <keys>. <keys> is a comma-separated list of keys. If None, the self._selected will be the focused one.
  --append <key>                                   Append Unimodal data <key> to any <keys> in --focus.

  --output-loom                                    Output loom-formatted file.
  --output-h5ad                                    Output h5ad-formatted file.

  --min-genes <number>                             Only keep cells with at least <number> of genes. [default: 500]
  --max-genes <number>                             Only keep cells with less than <number> of genes. [default: 6000]
  --min-umis <number>                              Only keep cells with at least <number> of UMIs.
  --max-umis <number>                              Only keep cells with less than <number> of UMIs.
  --mito-prefix <prefix>                           Prefix for mitochondrial genes. Can provide multiple prefixes for multiple organisms (e.g. "MT-" means to use "MT-", "GRCh38:MT-,mm10:mt-,MT-" means to use "MT-" for GRCh38, "mt-" for mm10 and "MT-" for all other organisms). [default: GRCh38:MT-,mm10:mt-,MT-]
  --percent-mito <percent>                         Only keep cells with mitochondrial percent less than <percent>%. [default: 20.0]
  --gene-percent-cells <percent>                   Only use genes that are expressed in at least <percent>% of cells to select variable genes. [default: 0.05]

  --output-filtration-results                      Output filtration results as a spreadsheet.
  --plot-filtration-results                        Plot filtration results as PDF files.
  --plot-filtration-figsize <figsize>              Figure size for filtration plots. <figsize> is a comma-separated list of two numbers, the width and height of the figure (e.g. 6,4).
  --min-genes-before-filtration <number>           If raw data matrix is input, empty barcodes will dominate pre-filtration statistics. To avoid this, for raw data matrix, only consider barcodes with at lease <number> genes for pre-filtration condition. [default: 100]

  --counts-per-cell-after <number>                 Total counts per cell after normalization. [default: 1e5]

  --select-hvf-flavor <flavor>                     Highly variable feature selection method. <flavor> can be 'pegasus' or 'Seurat'. [default: pegasus]
  --select-hvf-ngenes <nfeatures>                  Select top <nfeatures> highly variable features. If <flavor> is 'Seurat' and <nfeatures> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]
  --no-select-hvf                                  Do not select highly variable features.
  --plot-hvf                                       Plot highly variable feature selection.

  --pca-n <number>                                 Number of principal components. If Scanorama is used for batch correction, this parameter also sets Scanorama number of components. [default: 50]

  --random-state <seed>                            Random number generator seed. [default: 0]
  --temp-folder <temp_folder>                      Joblib temporary folder for memmapping numpy arrays.

  --calc-signature-scores <sig_list>               Calculate signature scores for gene sets in <sig_list>. <sig_list> is a comma-separated list of strings. Each string should either be a <GMT_file> or one of 'cell_cycle_human', 'cell_cycle_mouse', 'gender_human', 'gender_mouse', 'mitochondrial_genes_human', 'mitochondrial_genes_mouse', 'ribosomal_genes_human' and 'ribosomal_genes_mouse'.

  --correct-batch-effect                           Correct for batch effects.
  --batch <batch_attr>                             Use <batch_attr> as batches for batch correction. [default: Channel]
  --correction-method <method>                     Batch correction method, can be either 'harmony' for Harmony (Korsunsky et al. Nature Methods 2019), 'scanorama' for Scanorama (Hie et al. Nature Biotechnology 2019) or 'inmf' for integrative NMF (Yang and Michailidis Bioinformatics 2016, Welch et al. Cell 2019, Gao et al. Natuer Biotechnology 2021) [default: harmony]
  --harmony-nclusters <nclusters>                  Number of clusters used for Harmony batch correction.
  --inmf-lambda <lambda>                           Coefficient of regularization for iNMF. [default: 5.0]

  --nmf                                            Compute nonnegative matrix factorization (NMF) on highly variable features.
  --nmf-n <number>                                 Number of NMF components. IF iNMF is used for batch correction, this parameter also sets iNMF number of components. [default: 20]

  --knn-K <number>                                 Number of nearest neighbors for building kNN graph. [default: 100]
  --knn-full-speed                                 For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

  --kBET                                           Calculate kBET.
  --kBET-batch <batch>                             kBET batch keyword. [default: Channel]
  --kBET-alpha <alpha>                             kBET rejection alpha. [default: 0.05]
  --kBET-K <K>                                     kBET K. [default: 25]

  --diffmap                                        Calculate diffusion maps.
  --diffmap-ndc <number>                           Number of diffusion components. [default: 100]
  --diffmap-solver <solver>                        Solver for eigen decomposition, either 'eigsh' or 'randomized'. [default: eigsh]
  --diffmap-maxt <max_t>                           Maximum time stamp to search for the knee point. [default: 5000]
  --calculate-pseudotime <roots>                   Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  --louvain                                        Run louvain clustering algorithm.
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --louvain-class-label <label>                    Louvain cluster label name in result. [default: louvain_labels]

  --leiden                                         Run leiden clustering algorithm.
  --leiden-resolution <resolution>                 Resolution parameter for the leiden clustering algorithm. [default: 1.3]
  --leiden-niter <niter>                           Number of iterations of running the Leiden algorithm. If <niter> is negative, run Leiden iteratively until no improvement. [default: -1]
  --leiden-class-label <label>                     Leiden cluster label name in result. [default: leiden_labels]

  --spectral-louvain                               Run spectral-louvain clustering algorithm.
  --spectral-louvain-basis <basis>                 Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. If 'diffmap' is not calculated, use 'pca' instead. [default: diffmap]
  --spectral-louvain-nclusters <number>            Number of first level clusters for Kmeans. [default: 30]
  --spectral-louvain-nclusters2 <number>           Number of second level clusters for Kmeans. [default: 50]
  --spectral-louvain-ninit <number>                Number of Kmeans tries for first level clustering. Default is the same as scikit-learn Kmeans function. [default: 10]
  --spectral-louvain-resolution <resolution>       Resolution parameter for louvain. [default: 1.3]
  --spectral-louvain-class-label <label>           Spectral-louvain label name in result. [default: spectral_louvain_labels]

  --spectral-leiden                                Run spectral-leiden clustering algorithm.
  --spectral-leiden-basis <basis>                  Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. If 'diffmap' is not calculated, use 'pca' instead. [default: diffmap]
  --spectral-leiden-nclusters <number>             Number of first level clusters for Kmeans. [default: 30]
  --spectral-leiden-nclusters2 <number>            Number of second level clusters for Kmeans. [default: 50]
  --spectral-leiden-ninit <number>                 Number of Kmeans tries for first level clustering. Default is the same as scikit-learn Kmeans function. [default: 10]
  --spectral-leiden-resolution <resolution>        Resolution parameter for leiden. [default: 1.3]
  --spectral-leiden-class-label <label>            Spectral-leiden label name in result. [default: spectral_leiden_labels]

  --tsne                                           Run FIt-SNE package to compute t-SNE embeddings for visualization.
  --tsne-perplexity <perplexity>                   t-SNE perplexity parameter. [default: 30]
  --tsne-initialization <choice>                   <choice> can be either 'random' or 'pca'. 'random' refers to random initialization. 'pca' refers to PCA initialization as described in (CITE Kobak et al. 2019) [default: pca]

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

  --net-umap                                       Run net umap for visualization.
  --net-umap-polish-learning-rate <rate>           After running the deep regressor to predict new coordinate, what is the learning rate to use to polish the coordinates for UMAP. [default: 1.0]
  --net-umap-polish-nepochs <nepochs>              Number of iterations for polishing UMAP run. [default: 40]
  --net-umap-out-basis <basis>                     Output basis for net-UMAP. [default: net_umap]

  --net-fle                                        Run net FLE.
  --net-fle-polish-target-steps <steps>            After running the deep regressor to predict new coordinate, what is the number of force atlas 2 iterations. [default: 1500]
  --net-fle-out-basis <basis>                      Output basis for net-FLE. [default: net_fle]

  --infer-doublets                                 Infer doublets using the method described in https://github.com/klarman-cell-observatory/pegasus/raw/master/doublet_detection.pdf. Obs attribute 'doublet_score' stores Scrublet-like doublet scores and attribute 'demux_type' stores 'doublet/singlet' assignments.
  --expected-doublet-rate <rate>                   The expected doublet rate per sample. By default, calculate the expected rate based on number of cells from the 10x multiplet rate table.
  --dbl-cluster-attr <attr>                        <attr> refers to a cluster attribute containing cluster labels (e.g. 'louvain_labels'). Doublet clusters will be marked based on <attr> with the following criteria: passing the Fisher's exact test and having >= 50% of cells identified as doublets. By default, the first computed cluster attribute in the list of leiden, louvain, spectral_ledein and spectral_louvain is used.

  --citeseq                                        Input data contain both RNA and CITE-Seq modalities. This will set --focus to be the RNA modality and --append to be the CITE-Seq modality. In addition, 'ADT-' will be added in front of each antibody name to avoid name conflict with genes in the RNA modality.
  --citeseq-umap                                   For high quality cells kept in the RNA modality, generate a UMAP based on their antibody expression.
  --citeseq-umap-exclude <list>                    <list> is a comma-separated list of antibodies to be excluded from the UMAP calculation (e.g. Mouse-IgG1,Mouse-IgG2a).

  -h, --help                                       Print out help information.

Outputs:
  output_name.zarr.zip                     Output file in Zarr format. To load this file in python, use ``import pegasus; data = pegasus.read_input('output_name.zarr.zip')``. The log-normalized expression matrix is stored in ``data.X`` as a CSR-format sparse matrix. The ``obs`` field contains cell related attributes, including clustering results. For example, ``data.obs_names`` records cell barcodes; ``data.obs['Channel']`` records the channel each cell comes from; ``data.obs['n_genes']``, ``data.obs['n_counts']``, and ``data.obs['percent_mito']`` record the number of expressed genes, total UMI count, and mitochondrial rate for each cell respectively; ``data.obs['louvain_labels']`` and ``data.obs['approx_louvain_labels']`` record each cell's cluster labels using different clustring algorithms; ``data.obs['pseudo_time']`` records the inferred pseudotime for each cell. The ``var`` field contains gene related attributes. For example, ``data.var_names`` records gene symbols, ``data.var['gene_ids']`` records Ensembl gene IDs, and ``data.var['selected']`` records selected variable genes. The ``obsm`` field records embedding coordiates. For example, ``data.obsm['X_pca']`` records PCA coordinates, ``data.obsm['X_tsne']`` records tSNE coordinates, ``data.obsm['X_umap']`` records UMAP coordinates, ``data.obsm['X_diffmap']`` records diffusion map coordinates, and ``data.obsm['X_fle']`` records the force-directed layout coordinates from the diffusion components. The ``uns`` field stores other related information, such as reference genome (``data.uns['genome']``). This file can be loaded into R and converted into a Seurat object.
  output_name.<group>.h5ad                 Optional output. Only exists if '--output-h5ad' is set. Results in h5ad format per focused <group>. This file can be loaded into R and converted into a Seurat object.
  output_name.<group>.loom                 Optional output. Only exists if '--output-loom' is set. Results in loom format per focused <group>.
  output_name.<group>.filt.xlsx            Optional output. Only exists if '--output-filtration-results' is set. Filtration statistics per focused <group>. This file has two sheets --- Cell filtration stats and Gene filtration stats. The first sheet records cell filtering results and it has 10 columns: Channel, channel name; kept, number of cells kept; median_n_genes, median number of expressed genes in kept cells; median_n_umis, median number of UMIs in kept cells; median_percent_mito, median mitochondrial rate as UMIs between mitochondrial genes and all genes in kept cells; filt, number of cells filtered out; total, total number of cells before filtration, if the input contain all barcodes, this number is the cells left after '--min-genes-on-raw' filtration; median_n_genes_before, median expressed genes per cell before filtration; median_n_umis_before, median UMIs per cell before filtration; median_percent_mito_before, median mitochondrial rate per cell before filtration. The channels are sorted in ascending order with respect to the number of kept cells per channel. The second sheet records genes that failed to pass the filtering. This sheet has 3 columns: gene, gene name; n_cells, number of cells this gene is expressed; percent_cells, the fraction of cells this gene is expressed. Genes are ranked in ascending order according to number of cells the gene is expressed. Note that only genes not expressed in any cell are removed from the data. Other filtered genes are marked as non-robust and not used for TPM-like normalization.
  output_name.<group>.filt.gene.pdf        Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting gene count distributions before and after filtration per channel per focused <group>.
  output_name.<group>.filt.UMI.pdf         Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting UMI count distributions before and after filtration per channel per focused <group>.
  output_name.<group>.filt.mito.pdf        Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting mitochondrial rate distributions before and after filtration per channel per focused <group>.
  output_name.<group>.hvf.pdf              Optional output. Only exists if '--plot-hvf' is set. This file contains a scatter plot describing the highly variable gene selection procedure per focused <group>.
  output_name.<group>.<channel>.dbl.png    Optional output. Only exists if '--infer-doublets' is set. Each figure consists of 4 panels showing diagnostic plots for doublet inference. If there is only one channel in <group>, file name becomes output_name.<group>.dbl.png.

Examples:
  pegasus cluster -p 20 --correct-batch-effect --louvain --tsne manton_bm_10x.h5 manton_bm
  pegasus cluster -p 20 --leiden --umap --net-fle example.zarr.zip example_out
    """

    def execute(self):
        kwargs = {
            "n_jobs": int(self.args["--threads"]),
            "processed": self.args["--processed"],
            "black_list": self.args["--black-list"],
            "subcluster": False,
            "select_singlets": self.args["--select-singlets"],
            "remap_singlets": self.args["--remap-singlets"],
            "subset_singlets": self.args["--subset-singlets"],
            "genome": self.args["--genome"],
            "focus": self.split_string(self.args["--focus"]),
            "append": self.args["--append"],
            "output_h5ad": self.args["--output-h5ad"],
            "output_loom": self.args["--output-loom"],
            "min_genes": self.convert_to_int(self.args["--min-genes"]),
            "max_genes": self.convert_to_int(self.args["--max-genes"]),
            "min_umis": self.convert_to_int(self.args["--min-umis"]),
            "max_umis": self.convert_to_int(self.args["--max-umis"]),
            "mito_prefix": self.args["--mito-prefix"],
            "percent_mito": self.convert_to_float(self.args["--percent-mito"]),
            "percent_cells": float(self.args["--gene-percent-cells"]),
            "output_filt": self.args["<output_name>"]
            if self.args["--output-filtration-results"]
            else None,
            "plot_filt": self.args["<output_name>"]
            if self.args["--plot-filtration-results"]
            else None,
            "plot_filt_figsize": self.args["--plot-filtration-figsize"],
            "min_genes_before_filt": int(self.args["--min-genes-before-filtration"]),
            "norm_count": float(self.args["--counts-per-cell-after"]),
            "select_hvf": not self.args["--no-select-hvf"],
            "hvf_flavor": self.args["--select-hvf-flavor"],
            "hvf_ngenes": int(self.args["--select-hvf-ngenes"])
            if self.args["--select-hvf-ngenes"] != "None"
            else None,
            "plot_hvf": self.args["<output_name>"] if self.args["--plot-hvf"] else None,
            "batch_correction": self.args["--correct-batch-effect"],
            "batch_attr": self.args["--batch"],
            "correction_method": self.args["--correction-method"],
            "harmony_nclusters": self.convert_to_int(self.args["--harmony-nclusters"]),
            "inmf_lambda": self.convert_to_float(self.args["--inmf-lambda"]),
            "random_state": int(self.args["--random-state"]),
            "temp_folder": self.args["--temp-folder"],
            "calc_sigscore": self.args["--calc-signature-scores"],
            "pca_n": int(self.args["--pca-n"]),
            "nmf": self.args["--nmf"],
            "nmf_n": int(self.args["--nmf-n"]),
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
            "tsne_perplexity": float(self.args["--tsne-perplexity"]),
            "tsne_init": self.args["--tsne-initialization"],
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
            "infer_doublets": self.args["--infer-doublets"],
            "expected_doublet_rate": self.convert_to_float(self.args["--expected-doublet-rate"]),
            "dbl_cluster_attr": self.args["--dbl-cluster-attr"],
            "citeseq": self.args["--citeseq"],
            "citeseq_umap": self.args["--citeseq-umap"],
            "citeseq_umap_exclude": self.split_string(self.args["--citeseq-umap-exclude"]),
        }

        import logging
        logger = logging.getLogger("pegasus")
        fh = logging.FileHandler("{}.log".format(self.args["<output_name>"]))
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        run_pipeline(self.args["<input_file>"], self.args["<output_name>"], **kwargs)
