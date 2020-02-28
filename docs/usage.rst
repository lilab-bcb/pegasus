Use ``pegasus`` as a command line tool
---------------------------------------

``pegasus`` can be used as a command line tool. Type::

	pegasus -h

to see the help information::

	Usage:
		pegasus <command> [<args>...]
		pegasus -h | --help
		pegasus -v | --version

``pegasus`` has 14 sub-commands in 8 groups.

* Preprocessing:

	aggregate_matrix
		Aggregate sample count matrices into a single count matrix. It also enables users to import metadata into the count matrix.

* Demultiplexing:

	demuxEM
		Demultiplex cells/nuclei based on DNA barcodes for cell-hashing and nuclei-hashing data.

* Analyzing:
	
	cluster
		Perform first-pass analysis using the count matrix generated from 'aggregate_matrix'. This subcommand could perform low quality cell filtration, batch correction, variable gene selection, dimension reduction, diffusion map calculation, graph-based clustering, tSNE visualization. The final results will be written into h5ad-formatted file, which Seurat could load.
  		
	de_analysis
		Detect markers for each cluster by performing differential expression analysis per cluster (within cluster vs. outside cluster). DE tests include Welch's t-test, Fisher's exact test, Mann-Whitney U test. It can also calculate AUROC values for each gene.
    
	find_markers
		Find markers for each cluster by training classifiers using LightGBM.
    
	annotate_cluster
		This subcommand is used to automatically annotate cell types for each cluster based on existing markers. Currently, it works for human/mouse immune/brain cells.

* Plotting:

	plot
		Make static plots, which includes plotting tSNEs by cluster labels and different groups.
			
	iplot
		Make interactive plots using plotly. The outputs are HTML pages. You can visualize diffusion maps with this sub-command.

* Subclustering:

	view
		View attribute (e.g. cluster labels) and their values. This subcommand is used to determine cells to run subcluster analysis.

	subcluster
		Perform sub-cluster analyses on a subset of cells from the analyzed data (i.e. 'cluster' output).

* Web-based visualization:

	scp_output
		Generate output files for single cell portal.

	parquet
		Generate a PARQUET file for web-based visualization.	

* CITE-Seq:

	merge_rna_adt
		Merge RNA and ADT matrices into one 10x-formatted hdf5 file.

* MISC:

	check_indexes
		Check CITE-Seq/hashing indexes to avoid index collision.

---------------------------------


Quick guide
^^^^^^^^^^^

Suppose you have ``example.csv`` ready with the following contents::

	Sample,Source,Platform,Donor,Reference,Location
	sample_1,bone_marrow,NextSeq,1,GRCh38,/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5
	sample_2,bone_marrow,NextSeq,2,GRCh38,/my_dir/sample_2/filtered_gene_bc_matrices_h5.h5
	sample_3,pbmc,NextSeq,1,GRCh38,/my_dir/sample_3/filtered_gene_bc_matrices_h5.h5
	sample_4,pbmc,NextSeq,2,GRCh38,/my_dir/sample_4/filtered_gene_bc_matrices_h5.h5

You want to analyze all four samples but correct batch effects for bone marrow and pbmc samples separately. You can run the following commands::

	pegasus aggregate_matrix --attributes Source,Platform,Donor example.csv example
	pegasus cluster -p 20 --correct-batch-effect --batch-group-by Source -run-louvain --run-tsne example_10x.h5 example
	pegasus de_analysis --labels louvain_labels -p 20 --fisher example.h5ad example_de.xlsx
	pegasus annotate_cluster example.h5ad example.anno.txt
	pegasus plot composition --cluster-labels louvain_labels --attribute Donor --style normalized --not-stacked example.h5ad example.composition.pdf
	pegasus plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad example.scatter.pdf
	pegasus iplot --attribute louvain_labels diffmap_pca example.h5ad example.diffmap.html

The above analysis will give you tSNE, louvain cluster labels and diffusion maps in ``example.h5ad``. You can investigate donor-specific effects by looking at ``example.composition.pdf``. ``example.scatter.pdf`` plotted tSNE colored by louvain_labels and Donor info side-by-side. You can explore the diffusion map in 3D by looking at ``example.diffmap.html``. This html maps all diffusion components into 3D using PCA.

If you want to perform subcluster analysis by combining cluster 1 and 3, run the following command::

	pegasus subcluster -p 20 --correct-batch-effect example.h5ad 1,3 example_sub


---------------------------------


``pegasus aggregate_matrix``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step for single cell analysis is to generate one count matrix from cellranger's channel-specific count matrices. ``pegasus aggregate_matrix`` allows aggregating arbitrary matrices with the help of a *CSV* file.

Type::

	pegasus aggregate_matrix -h

to see the usage information::

	Usage:
		pegasus aggregate_matrix <csv_file> <output_name> [--restriction <restriction>... --attributes <attributes> --default-reference <reference> --select-only-singlets --minimum-number-of-genes <ngene>]
		pegasus aggregate_matrix -h

* Arguments:

	csv_file
		Input csv-formatted file containing information of each sc/snRNA-seq sample. This file must contain at least 2 columns - Sample, sample name and Location, location of the sample count matrix in either 10x v2/v3, DGE, mtx, csv, tsv or loom format. Additionally, an optional Reference column can be used to select samples generated from a same reference (e.g. mm10). If the count matrix is in either DGE, mtx, csv, tsv, or loom format, the value in this column will be used as the reference since the count matrix file does not contain reference name information. In addition, the Reference column can be used to aggregate count matrices generated from different genome versions or gene annotations together under a unified reference. For example, if we have one matrix generated from mm9 and the other one generated from mm10, we can write mm9_10 for these two matrices in their Reference column. Pegasus will change their references to 'mm9_10' and use the union of gene symbols from the two matrices as the gene symbols of the aggregated matrix. For HDF5 files (e.g. 10x v2/v3), the reference name contained in the file does not need to match the value in this column. In fact, we use this column to rename references in HDF5 files. For example, if we have two HDF files, one generated from mm9 and the other generated from mm10. We can set these two files' Reference column value to 'mm9_10', which will rename their reference names into mm9_10 and the aggregated matrix will contain all genes from either mm9 or mm10. This renaming feature does not work if one HDF5 file contain multiple references (e.g. mm10 and GRCh38). See below for an example csv::

			Sample,Source,Platform,Donor,Reference,Location
 			sample_1,bone_marrow,NextSeq,1,GRCh38,/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5
			sample_2,bone_marrow,NextSeq,2,GRCh38,/my_dir/sample_2/filtered_gene_bc_matrices_h5.h5
			sample_3,pbmc,NextSeq,1,GRCh38,/my_dir/sample_3/filtered_gene_bc_matrices_h5.h5
			sample_4,pbmc,NextSeq,2,GRCh38,/my_dir/sample_4/filtered_gene_bc_matrices_h5.h5

	output_name
		The output file name.

* Options:
	
	-\\-restriction <restriction>...
		Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value or name:~value,..,value, where ~ refers to not. You can specifiy multiple restrictions by setting this option multiple times.

	-\\-attributes <attributes>
		Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.

	-\\-default-reference <reference>
		If sample count matrix is in either DGE, mtx, csv, tsv or loom format and there is no Reference column in the csv_file, use <reference> as the reference.

	-\\-select-only-singlets
		If we have demultiplexed data, turning on this option will make pegasus only include barcodes that are predicted as singlets.

	-\\-minimum-number-of-genes <ngene>
		Only keep barcodes with at least <ngene> expressed genes.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.h5sc
		A pegasus-formatted HDF5 file containing the count matrices and associated attributes.

* Examples::

	pegasus aggregate_matrix --restriction Source:BM,CB --restriction Individual:1-8 --attributes Source,Platform Manton_count_matrix.csv manton_bm_cb


---------------------------------


``pegasus demuxEM``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have data generated by cell-hashing or nuclei-hashing, you can use ``pegasus demuxEM`` to demultiplex your data. 

Type::

	pegasus demuxEM -h

to see the usage information::

	Usage:
		pegasus demuxEM [options] <input_adt_csv_file> <input_raw_gene_bc_matrices_h5.h5> <output_name>
		pegasus demuxEM -h

* Arguments:

	input_adt_csv_file
		Input ADT (antibody tag) count matrix in CSV format.

	input_raw_gene_bc_matrices_h5.h5
		Input raw RNA expression matrix in 10x hdf5 format.

	output_name
		Output name. All outputs will use it as the prefix.

* Options:

  	\-p <number>, -\\-threads <number>
		Number of threads. [default: 1]

	-\\-genome <genome>
		Reference genome name. If not provided, we will infer it from the expression matrix file.

	-\\-alpha-on-samples <alpha>
		The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse. [default: 0.0]

	-\\-min-num-genes <number>
		We only demultiplex cells/nuclei with at least <number> of expressed genes. [default: 100]

	-\\-min-num-umis <number>
		We only demultiplex cells/nuclei with at least <number> of UMIs. [default: 100] 

	-\\-min-signal-hashtag <count>
		Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown. [default: 10.0]

	-\\-random-state <seed>
		The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]

	-\\-generate-diagnostic-plots
		Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc.

	-\\-generate-gender-plot <genes>
		Generate violin plots using gender-specific genes (e.g. Xist). <gene> is a comma-separated list of gene names. 
	
	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name_demux.h5sc
		RNA expression matrix with demultiplexed sample identities in pegasus HDF5 format.

	output_name_ADTs.h5ad
		Antibody tag matrix in h5ad format.

	output_name_demux.h5ad
		Demultiplexed RNA count matrix in h5ad format.

	output_name.ambient_hashtag.hist.pdf
		Optional output. A histogram plot depicting hashtag distributions of empty droplets and non-empty droplets.

	output_name.background_probabilities.bar.pdf
		Optional output. A bar plot visualizing the estimated hashtag background probability distribution.

	output_name.real_content.hist.pdf
		Optional output. A histogram plot depicting hashtag distributions of not-real-cells and real-cells as defined by total number of expressed genes in the RNA assay.

	output_name.rna_demux.hist.pdf
		Optional output. A histogram plot depicting RNA UMI distribution for singlets, doublets and unknown cells.

	output_name.gene_name.violin.pdf
		Optional outputs. Violin plots depicting gender-specific gene expression across samples. We can have multiple plots if a gene list is provided in '--generate-gender-plot' option.

* Examples::

	pegasus demuxEM -p 8 --hash-type cell-hashing --generate-diagnostic-plots example_adt.csv example_raw_gene_bc_matrices_h5.h5 example_output


---------------------------------


``pegasus cluster``
^^^^^^^^^^^^^^^^^^^

Once we collected the count matrix in 10x (``example_10x.h5``) or pegasus (``example.h5sc``) format, we can perform single cell analysis using ``pegasus cluster``.

Type::

	pegasus cluster -h

to see the usage information::

	Usage:
		pegasus cluster [options] <input_file> <output_name>
		pegasus cluster -h

* Arguments:

	input_file
		Input file in 10x or pegasus format. If first-pass analysis has been performed, but you want to run some additional analysis, you could also pass a h5ad-formatted file.

	output_name      
		Output file name. All outputs will use it as the prefix.

* Options:

	\-p <number>, -\\-threads <number>
		Number of threads. [default: 1]

	-\\-processed
		Input file is processed and thus no PCA & diffmap will be run.

	-\\-considered-refs <ref_list>
		A string contains comma-separated reference(e.g. genome) names. pegasus will read all groups associated with reference names in the list from the input file. If <ref_list> is None, all groups will be considered. For formats like loom, mtx, dge, csv and tsv, genome is used to provide genome name. In this case if genome is None, except mtx format, '' is used as genome name instead.
  
  	-\\-channel <channel_attr>
		Use <channel_attr> to create a 'Channel' column metadata field. All cells within a channel are assumed to come from a same batch.

	-\\-black-list <black_list>
		Cell barcode attributes in black list will be popped out. Format is "attr1,attr2,...,attrn".

	-\\-min-genes-on-raw <number>
		If input are raw 10x matrix, which include all barcodes, perform a pre-filtration step to keep the data size small. In the pre-filtration step, only keep cells with at least <number> of genes. [default: 100]

	-\\-select-singlets
		Only select DemuxEM-predicted singlets for analysis.

	-\\-remap-singlets <remap_string> 
	Remap singlet names using <remap_string>, where <remap_string> takes the format "new_name_i:old_name_1,old_name_2;new_name_ii:old_name_3;...". For example, if we hashed 5 libraries from 3 samples sample1_lib1, sample1_lib2, sample2_lib1, sample2_lib2 and sample3, we can remap them to 3 samples using this string: "sample1:sample1_lib1,sample1_lib2;sample2:sample2_lib1,sample2_lib2". After that, original singlet names will be kept in metadate field with key name 'assignment.orig'.

	-\\-subset-singlets <subset_string>
	If select singlets, only select singlets in the <subset_string>, which takes the format "name1,name2,...". Note that if --remap-singlets is specified, subsetting happens after remapping. For example, we can only select singlets from sampe 1 and 3 using "sample1,sample3".

	-\\-cite-seq
		Data are CITE-Seq data. pegasus will perform analyses on RNA count matrix first. Then it will attach the ADT matrix to the RNA matrix with all antibody names changing to 'AD-' + antibody_name. Lastly, it will embed the antibody expression using FIt-SNE (the basis used for plotting is 'citeseq_fitsne').

	-\\-cite-seq-capping <percentile>
		For CITE-Seq surface protein expression, make all cells with expression > <percentile> to the value at <percentile> to smooth outlier. Set <percentile> to 100.0 to turn this option off. [default: 99.99]

  	-\\-output-filtration-results
		Output filtration results as a spreadsheet.

	-\\-plot-filtration-results
		Plot filtration results as PDF files.

	-\\-plot-filtration-figsize <figsize>
		Figure size for filtration plots. <figsize> is a comma-separated list of two numbers, the width and height of the figure (e.g. 6,4).

	-\\-output-seurat-compatible
		Output seurat-compatible h5ad file. Caution: File size might be large, do not turn this option on for large data sets.

	-\\-output-loom
		Output loom-formatted file.

  	-\\-min-genes <number>
		Only keep cells with at least <number> of genes. [default: 500]

	-\\-max-genes <number>
		Only keep cells with less than <number> of genes. [default: 6000]

	-\\-min-umis <number>
		Only keep cells with at least <number> of UMIs. [default: 100]

	-\\-max-umis <number>
		Only keep cells with less than <number> of UMIs. [default: 600000]

	-\\-mito-prefix <prefix>
		Prefix for mitochondrial genes. If multiple prefixes are provided, separate them by comma (e.g. "MT-,mt-"). [default: MT-]

	-\\-percent-mito <ratio>
		Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]

	-\\-gene-percent-cells <ratio>
		Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]

	-\\-counts-per-cell-after <number>
		Total counts per cell after normalization. [default: 1e5]

	-\\-select-hvf-flavor <flavor>
		Highly variable feature selection method. <flavor> can be 'pegasus' or 'Seurat'. [default: pegasus]

	-\\-select-hvf-ngenes <nfeatures>
		Select top <nfeatures> highly variable features. If <flavor> is 'Seurat' and <ngenes> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]

	-\\-no-select-hvf
		Do not select highly variable features.

	-\\-plot-hvf
		Plot highly variable feature selection.

	-\\-correct-batch-effect
		Correct for batch effects.

	-\\-correction-method <method>
		Batch correction method, can be either 'L/S' for location/scale adjustment algorithm (Li and Wong. The analysis of Gene Expression Data 2003) or 'harmony' for Harmony (Korsunsky et al. Nature Methods 2019). [default: harmony]

	-\\-batch-group-by <expression>
		Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+...+attrn', or 'attr=value11,...,value1n_1;value21,...,value2n_2;...;valuem1,...,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',...,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,...,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.

	-\\-random-state <seed>
		Random number generator seed. [default: 0]

	-\\-temp-folder <temp_folder>
		Joblib temporary folder for memmapping numpy arrays.
  
	-\\-nPC <number>
		Number of principal components. [default: 50]

	-\\-knn-K <number>
		Number of nearest neighbors for building kNN graph. [default: 100]

	-\\-knn-full-speed
		For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

	-\\-kBET
		Calculate kBET.

	-\\-kBET-batch <batch>
		kBET batch keyword.

	-\\-kBET-alpha <alpha>
		kBET rejection alpha. [default: 0.05]

	-\\-kBET-K <K>
		kBET K. [default: 25]

	-\\-diffmap
		Calculate diffusion maps.

	-\\-diffmap-ndc <number>
		Number of diffusion components. [default: 50]

	-\\-diffmap-alpha <alpha>
		Power parameter for diffusion-based pseudotime. [default: 0.5]

	-\\-diffmap-solver <solver>
		Solver for eigen decomposition, either 'randomized' or 'eigsh'. [default: randomized]

	-\\-diffmap-to-3d
		If map diffusion map into 3D space using PCA.

	-\\-calculate-pseudotime <roots>
		Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  	-\\-louvain
  		Run louvain clustering algorithm.

	-\\-louvain-resolution <resolution>
		Resolution parameter for the louvain clustering algorithm. [default: 1.3]

	-\\-louvain-class-label <label>
		Louvain cluster label name in AnnData. [default: louvain_labels]

	-\\-leiden
		Run leiden clustering algorithm.

	-\\-leiden-resolution <resolution>
		Resolution parameter for the leiden clustering algorithm. [default: 1.3]

	-\\-leiden-niter <niter>
		Number of iterations of running the Leiden algorithm. If <niter> is negative, run Leiden iteratively until no improvement. [default: -1]

	-\\-leiden-class-label <label>
		Leiden cluster label name in AnnData. [default: leiden_labels]

	-\\-spectral-louvain
		Run spectral-louvain clustering algorithm.

	-\\-spectral-louvain-basis <basis>
		Basis used for KMeans clustering. Can be 'pca', or 'diffmap'. [default: diffmap]

	-\\-spectral-louvain-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-spectral-louvain-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-spectral-louvain-resolution <resolution>.
		Resolution parameter for louvain. [default: 1.3]

	-\\-spectral-louvain-class-label <label>
		Spectral-louvain label name in AnnData. [default: spectral_louvain_labels]

	-\\-spectral-leiden
		Run spectral-leiden clustering algorithm.

	-\\-spectral-leiden-basis <basis>
		Basis used for KMeans clustering. Can be 'pca', or 'diffmap'. [default: diffmap]

	-\\-spectral-leiden-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-spectral-leiden-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-spectral-leiden-resolution <resolution>
		Resolution parameter for leiden. [default: 1.3]

	-\\-spectral-leiden-class-label <label>
		Spectral-leiden label name in AnnData. [default: spectral_leiden_labels]

	-\\-tsne
		Run multi-core t-SNE for visualization.

	-\\-fitsne
  		Run FIt-SNE for visualization.

	-\\-tsne-perplexity <perplexity>
		t-SNE's perplexity parameter, used by both tSNE, FItSNE and net-tSNE. [default: 30]

  	-\\-umap
  		Run umap for visualization.

	-\\-umap-K <K>
		K neighbors for umap. [default: 15]

	-\\-umap-min-dist <number>
		Umap parameter. [default: 0.5]

	-\\-umap-spread <spread>
		Umap parameter. [default: 1.0]

	-\\-fle
		Run force-directed layout embedding.

	-\\-fle-K <K>
		K neighbors for building graph for FLE. [default: 50]

	-\\-fle-target-change-per-node <change>
		Target change per node to stop forceAtlas2. [default: 2.0]

	-\\-fle-target-steps <steps>
		Maximum number of iterations before stopping the forceAtlas2 algoritm. [default: 5000]

	-\\-fle-memory <memory>
		Memory size in GB for the Java FA2 component. [default: 8]

	-\\-net-down-sample-fraction <frac>
		Down sampling fraction for net-related visualization. [default: 0.1]

	-\\-net-down-sample-K <K>
		Use <K> neighbors to estimate local density for each data point for down sampling. [default: 25]

	-\\-net-down-sample-alpha <alpha>
		Weighted down sample, proportional to radius^alpha. [default: 1.0]

	-\\-net-regressor-L2-penalty <value>
		L2 penalty parameter for the deep net regressor. [default: 0.1]

	-\\-net-tsne
		Run net tSNE for visualization.

	-\\-net-tsne-polish-learning-frac <frac>
		After running the deep regressor to predict new coordinates, use <frac> * nsample as the learning rate to use to polish the coordinates. [default: 0.33]

	-\\-net-tsne-polish-niter <niter>
		Number of iterations for polishing tSNE run. [default: 150]

	-\\-net-tsne-out-basis <basis>
		Output basis for net-tSNE. [default: net_tsne]

	-\\-run-net-umap
		Run net umap for visualization.

	-\\-net-umap-polish-learning-rate <rate>
		After running the deep regressor to predict new coordinate, what is the learning rate to use to polish the coordinates for UMAP. [default: 1.0]

	-\\-net-umap-polish-nepochs <nepochs>
		Number of iterations for polishing UMAP run. [default: 40]

	-\\-net-umap-out-basis <basis>
		Output basis for net-UMAP. [default: net_umap]

	-\\-net-fle
		Run net FLE.

	-\\-net-fle-polish-target-steps <steps>
		After running the deep regressor to predict new coordinate, what is the number of force atlas 2 iterations. [default: 1500]

	-\\-net-fle-out-basis <basis>
		Output basis for net-FLE. [default: net_fle]

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.h5ad
		Output file in h5ad format. To load this file in python, use ``import pegasus; data = pegasus.tools.read_input('output_name.h5ad', mode = 'a')``. The log-normalized expression matrix is stored in ``data.X`` as a CSR-format sparse matrix. The ``obs`` field contains cell related attributes, including clustering results. For example, ``data.obs_names`` records cell barcodes; ``data.obs['Channel']`` records the channel each cell comes from; ``data.obs['n_genes']``, ``data.obs['n_counts']``, and ``data.obs['percent_mito']`` record the number of expressed genes, total UMI count, and mitochondrial rate for each cell respectively; ``data.obs['louvain_labels']`` and ``data.obs['approx_louvain_labels']`` record each cell's cluster labels using different clustring algorithms; ``data.obs['pseudo_time']`` records the inferred pseudotime for each cell. The ``var`` field contains gene related attributes. For example, ``data.var_names`` records gene symbols, ``data.var['gene_ids']`` records Ensembl gene IDs, and ``data.var['selected']`` records selected variable genes. The ``obsm`` field records embedding coordiates. For example, ``data.obsm['X_pca']`` records PCA coordinates, ``data.obsm['X_tsne']`` records tSNE coordinates, ``data.obsm['X_umap']`` records UMAP coordinates, ``data.obsm['X_diffmap']`` records diffusion map coordinates, ``data.obsm['X_diffmap_pca']`` records the first 3 PCs by projecting the diffusion components using PCA, and ``data.obsm['X_fle']`` records the force-directed layout coordinates from the diffusion components. The ``uns`` field stores other related information, such as reference genome (``data.uns['genome']``). If '--make-output-seurat-compatible' is on, this file can be loaded into R and converted into a Seurat object.

	output_name.seruat.h5ad
		Optional output. Only exists if '--output-seruat-compatible' is set. 'output_name.h5ad' in seurat-compatible manner. This file can be loaded into R and converted into a Seurat object.

	output_name.loom
		Optional output. Only exists if '--output-loom' is set. 'output_name.h5ad' in loom format for visualization.

	output_name.filt.xlsx
		Optional output. Only exists if '--output-filtration-results' is set. This file has two sheets --- Cell filtration stats and Gene filtration stats. The first sheet records cell filtering results and it has 10 columns: Channel, channel name; kept, number of cells kept; median_n_genes, median number of expressed genes in kept cells; median_n_umis, median number of UMIs in kept cells; median_percent_mito, median mitochondrial rate as UMIs between mitochondrial genes and all genes in kept cells; filt, number of cells filtered out; total, total number of cells before filtration, if the input contain all barcodes, this number is the cells left after '--min-genes-on-raw' filtration; median_n_genes_before, median expressed genes per cell before filtration; median_n_umis_before, median UMIs per cell before filtration; median_percent_mito_before, median mitochondrial rate per cell before filtration. The channels are sorted in ascending order with respect to the number of kept cells per channel. The second sheet records genes that failed to pass the filtering. This sheet has 3 columns: gene, gene name; n_cells, number of cells this gene is expressed; percent_cells, the fraction of cells this gene is expressed. Genes are ranked in ascending order according to number of cells the gene is expressed. Note that only genes not expressed in any cell are removed from the data. Other filtered genes are marked as non-robust and not used for TPM-like normalization.

	output_name.filt.gene.pdf
		Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting gene count distributions before and after filtration per channel.
	
	output_name.filt.UMI.pdf
		Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting UMI count distributions before and after filtration per channel.
	
	output_name.filt.mito.pdf
		Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting mitochondrial rate distributions before and after filtration per channel.

* Examples::

	pegasus cluster -p 20 --correct-batch-effect --louvain --tsne example_10x.h5 example
	pegasus cluster -p 20 --leiden --umap --net-fle example.h5sc example


---------------------------------


``pegasus de_analysis``
^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the clusters, we can detect markers using ``pegasus de_analysis``.

Type::

	pegasus de_analysis -h

to see the usage information::

	Usage:
		pegasus de_analysis [options] <input_h5ad_file> <output_spreadsheet>
		pegasus de_analysis -h

* Arguments:

	input_h5ad_file
		Single cell data with clustering calculated. DE results would be written back.
	
	output_spreadsheet
		Output spreadsheet with DE results.

* Options:

	\-p <threads>
		Use <threads> threads. [default: 1]

	-\\-labels <attr>
		<attr> used as cluster labels. [default: louvain_labels]

	-\\-result-key <key>
		Store DE results into AnnData varm with key = <key>. [default: de_res]

	-\\-auc
		Calculate area under ROC (AUROC) and area under Precision-Recall (AUPR).

	-\\-t
		Calculate Welch's t-test.

	-\\-fisher
		Calculate Fisher's exact test.

	-\\-mwu
		Calculate Mann-Whitney U test.

	-\\-temp-folder <temp_folder>
		Joblib temporary folder for memmapping numpy arrays.

	-\\-alpha <alpha>
		Control false discovery rate at <alpha>. [default: 0.05]

	-\\-ndigits <ndigits>
		Round non p-values and q-values to <ndigits> after decimal point in the excel. [default: 3]

	-\\-quiet 
		Do not show detailed intermediate outputs.

	\-h, -\\-help
		Print out help information.

* Outputs:

	input_h5ad_file
		DE results would be written back to the 'varm' field with name set by '--result-key <key>'.

	output_spreadsheet
		An excel spreadsheet containing DE results. Each cluster has two tabs in the spreadsheet. One is for up-regulated genes and the other is for down-regulated genes.

* Examples::

	pegasus de_analysis -p 26 --labels louvain_labels --auc --t --fisher --mwu example.h5ad example_de.xlsx


---------------------------------


``pegasus find_markers``
^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the DE results, we can optionally find cluster-specific markers with gradient boosting using ``pegasus find_markers``.

Type::

	pegasus find_markers -h

to see the usage information::

	Usage:
		pegasus find_markers [options] <input_h5ad_file> <output_spreadsheet>
		pegasus find_markers -h

* Arguments:

	input_h5ad_file
		Single cell data after running the de_analysis.

	output_spreadsheet
		Output spreadsheet with LightGBM detected markers.

* Options:

	\-p <threads>
		Use <threads> threads. [default: 1]

	-\\-labels <attr>
		<attr> used as cluster labels. [default: louvain_labels]

	-\\-de_key <key>
		Key for storing DE results in 'varm' field.

	-\\-remove-ribo
		Remove ribosomal genes with either RPL or RPS as prefixes.

	-\\-min-gain <gain>
		Only report genes with a feature importance score (in gain) of at least <gain>. [default: 1.0]

	-\\-random-state <seed>
		Random state for initializing LightGBM and KMeans. [default: 0]

	

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_spreadsheet
		An excel spreadsheet containing detected markers. Each cluster has one tab in the spreadsheet and each tab has six columns, listing markers that are strongly up-regulated, weakly up-regulated, down-regulated and their associated LightGBM gains.

* Examples::

	pegasus find_markers --labels louvain_labels --remove-ribo --min-gain 10.0 -p 10 example.h5ad example.markers.xlsx


---------------------------------


``pegasus annotate_cluster``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the DE results, we could optionally identify putative cell types for each cluster using ``pegasus annotate_cluster``. Currently, this subcommand works for human/mouse immune/brain cells. This command has two forms: the first form generates putative annotations and the second form write annotations into the h5ad object.

Type::

	pegasus annotate_cluster -h

to see the usage information::

	Usage:
		pegasus annotate_cluster [--marker-file <file> --de-test <test> --de-alpha <alpha> --de-key <key> --minimum-report-score <score> --do-not-use-non-de-genes] <input_h5ad_file> <output_file>
		pegasus annotate_cluster --annotation <annotation_string> <input_h5ad_file>
		pegasus annotate_cluster -h

* Arguments:

	input_h5ad_file
		Single cell data with DE analysis done by ``pegasus de_analysis``.

	output_file
		Output annotation file.

* Options:

	-\\-marker-file <file>
		JSON file for markers. Could also be ``human_immune``/``mouse_immune``/``mouse_brain``/``human_brain``, which triggers pegasus to markers included in the package. [default: human_immune]

	-\\-de-test <test>
		DE test to use to infer cell types. [default: t]

	-\\-de-alpha <alpha>
		False discovery rate to control family-wise error rate. [default: 0.05]

	-\\-de-key <key>
		Keyword where the DE results store in 'varm' field. [default: de_res]

	-\\-minimum-report-score <score>
		Minimum cell type score to report a potential cell type. [default: 0.5]

	-\\-do-not-use-non-de-genes
		Do not count non DE genes as down-regulated.

	-\\-annotation <annotation_string>
		Write cell type annotations in <annotation_string> into <input_h5ad_file>. <annotation_string> has this format: ``'anno_name:clust_name:anno_1;anno_2;...;anno_n'``, where ``anno_name`` is the annotation attribute in the h5ad object, ``clust_name`` is the attribute with cluster ids, and ``anno_i`` is the annotation for cluster i.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_file
		This is a text file. For each cluster, all its putative cell types are listed in descending order of the cell type score. For each putative cell type, all markers support this cell type are listed. If one putative cell type has cell subtypes, all subtypes will be listed under this cell type.

* Examples::

	pegasus annotate_cluster example.h5ad example.anno.txt
	pegasus annotate_cluster --annotation "anno:louvain_labels:T cells;B cells;NK cells;Monocytes" example.h5ad


---------------------------------



``pegasus plot``
^^^^^^^^^^^^^^^^^

We can make a variety of figures using ``pegasus plot``.

Type::

	pegasus plot -h

to see the usage information::

	Usage:
  		pegasus plot [options] [--restriction <restriction>...] <plot_type> <input_h5ad_file> <output_file>
		pegasus plot -h

* Arguments:

	plot_type
		Only 2D plots, chosen from 'composition', 'scatter', 'scatter_groups', 'scatter_genes', 'scatter_gene_groups', 'heatmap', and 'qc_violin'.

	input_h5ad_file
		Single cell data in h5ad file format with clustering done by ``pegasus cluster``.

  	output_file
  		Output image file.

* Options:

	-\\-dpi <dpi>
		DPI value for the figure. [default: 500]

	-\\-cluster-labels <attr>
		Use <attr> as cluster labels. This option is used in 'composition', 'scatter_groups', 'heatmap', and 'qc_violin'.

  	-\\-attribute <attr>
  		Plot <attr> against cluster labels. This option is only used in 'composition' and 'qc_violin'.

	-\\-basis <basis>
		Basis for 2D plotting, chosen from 'tsne', 'fitsne', 'umap', 'pca', 'rpca', 'fle', 'diffmap_pca', 'net_tsne', 'net_fitsne', 'net_umap' or 'net_fle'. If CITE-Seq data is used, basis can also be 'citeseq_fitsne'. This option is used in 'scatter', 'scatter_groups', 'scatter_genes', and 'scatter_gene_groups'. [default: fitsne]

	-\\-attributes <attrs>
		<attrs> is a comma-separated list of attributes to color the basis. This option is only used in 'scatter'.

	-\\-restriction <restriction>...
		Set restriction if you only want to plot a subset of data. Multiple <restriction> strings are allowed. Each <restriction> takes the format of 'attr:value,value', or 'attr:~value,value..' which means excluding values. This option is used in 'composition' and 'scatter'.
	
	-\\-apply-to-each-figure
		Indicate that the <restriction> strings are not applied to all attributes but for specific attributes. The string's 'attr' value should math the attribute you want to restrict. 

	-\\-show-background
		Show points that are not selected as gray.

	-\\-group <attr>
		<attr> is used to make group plots. In group plots, the first one contains all components in the group and the following plots show each component separately. This option is iused in 'scatter_groups' and 'scatter_gene_groups'. If <attr> is a semi-colon-separated string, parse the string as groups.

	-\\-genes <genes>
		<genes> is a comma-separated list of gene names to visualize. This option is used in 'scatter_genes' and 'heatmap'.

	-\\-gene <gene>
  		Visualize <gene> in group plots. This option is only used in 'scatter_gene_groups'.

	-\\-style <style>
		Composition plot styles. Can be either 'frequency', 'count', or 'normalized'. [default: frequency]

	-\\-not-stacked
		Do not stack bars in composition plot.
  
	-\\-log-y
		Plot y axis in log10 scale for composition plot.

	-\\-nrows <nrows>
		Number of rows in the figure. If not set, pegasus will figure it out automatically.

	-\\-ncols <ncols>
		Number of columns in the figure. If not set, pegasus will figure it out automatically.

	-\\-subplot-size <sizes>
		Sub-plot size in inches, w x h, separated by comma. Note that margins are not counted in the sizes. For composition, default is (6, 4). For scatter plots, default is (4, 4).

	-\\-left <left>
		Figure's left margin in fraction with respect to subplot width.

	-\\-bottom <bottom>
		Figure's bottom margin in fraction with respect to subplot height.

	-\\-wspace <wspace>
		Horizontal space between subplots in fraction with respect to subplot width.

	-\\-hspace <hspace>
		Vertical space between subplots in fraction with respect to subplot height.

	-\\-alpha <alpha>
		Point transparent parameter.

	-\\-legend-fontsize <fontsize>
		Legend font size.

	-\\-use-raw
		Use anndata stored raw expression matrix. Only used by 'scatter_genes' and 'scatter_gene_groups'.

	-\\-do-not-show-all
		Do not show all components in group for scatter_groups.

	-\\-show-zscore
		If show zscore in heatmap.

	-\\-heatmap-title <title>
		Title for heatmap.

	-\\-qc-type <type>
		Plot qc_violin by annotation, <type> can be either 'gene', 'count' (UMI), or 'mito' (mitochondrial rate). [default: gene]

	-\\-qc-xtick-font <font>
		x tick font for qc_violin. [default: 5]

	-\\-qc-xtick-rotation
		If rorate x label.

	-\\-qc-line-width <width>
		Line width for qc_violin. [default: 0.5]

	\-h, -\\-help
		Print out help information.

Examples::

	pegasus plot composition --cluster-labels louvain_labels --attribute Donor --style normalized --not-stacked example.h5ad example.composition.pdf
	pegasus plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad example.scatter.pdf
	pegasus plot scatter_groups --cluster-labels louvain_labels --group Donor example.h5ad example.scatter_groups.pdf
	pegasus plot scatter_genes --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP example.h5ad example.genes.pdf
	pegasus plot scatter_gene_groups --gene CD8A --group Donor example.h5ad example.gene_groups.pdf
	pegasus plot heatmap --cluster-labels louvain_labels --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP --heatmap-title 'markers' example.h5ad example.heatmap.pdf
	pegasus plot qc_violin --qc-type gene --cluster-labels louvain_labels --attribute Channel --subplot-size 7,5 --qc-xtick-font 5 --qc-line-width 0.5 example.h5ad example.qc_violin.pdf


---------------------------------


``pegasus iplot``
^^^^^^^^^^^^^^^^^^

We can also make interactive plots in html format using ``pegasus iplot``. These interactive plots are very helpful if you want to explore the diffusion maps.

Type::

	pegasus iplot -h

to see the usage information::

	Usage:
		pegasus iplot --attribute <attr> [options] <basis> <input_h5ad_file> <output_html_file>
		pegasus iplot -h

* Arguments:

	basis
		Basis can be either 'tsne', 'fitsne', 'umap', 'diffmap', 'pca', or 'diffmap_pca'.
	
	input_h5ad_file
		Single cell data with clustering done in h5ad file format.

	output_html_file
		Output interactive plot in html format.

* Options:

	-\\-attribute <attr>
		Use attribute <attr> as labels in the plot.

	-\\-is-real
		<attr> is real valued.

	-\\-is-gene
		<attr> is a gene name.

	-\\-log10
		If take log10 of real values.

	\-h, -\\-help
		Print out help information.

* Examples::

	pegasus iplot --attribute louvain_labels tsne example.h5ad example.tsne.html
	pegasus iplot --attribute louvain_labels diffmap_pca example.h5ad example.diffmap.html


---------------------------------


``pegasus view``
^^^^^^^^^^^^^^^^^

We may want to further perform sub-cluster analysis on a subset of cells. This sub-command helps us to define the subset.

Type::

	pegasus view -h

to see the usage information::

	Usage:
		pegasus view [--show-attributes --show-gene-attributes --show-values-for-attributes <attributes>] <input_h5ad_file>
		pegasus view -h

* Arguments:

	input_h5ad_file
		Analyzed single cell data in h5ad format.

* Options:

	-\\-show-attributes
  		Show the available sample attributes in the input dataset.

	-\\-show-gene-attributes
		Show the available gene attributes in the input dataset.

	-\\-show-values-for-attributes <attributes>
		Show the available values for specified attributes in the input dataset. <attributes> should be a comma-separated list of attributes.

	\-h, -\\-help
		Print out help information.

* Examples::

	pegasus view --show-attributes example.h5ad
	pegasus view --show-gene-attributes example.h5ad
	pegasus view --show-values-for-attributes louvain_labels,Donor example.h5ad


---------------------------------


``pegasus subcluster``
^^^^^^^^^^^^^^^^^^^^^^^

If there is a subset of cells that we want to further cluster, we can run ``pegasus subcluster``. This sub-command will outputs a new h5ad file that you can run ``de_analysis``, ``plot`` and ``iplot`` on.

Type::

	pegasus subcluster -h

to see the usage information::

	Usage:
		pegasus subcluster [options] --subset-selection <subset-selection>... <input_file> <output_name>
		pegasus subcluster -h

* Arguments:

	input_file
		Single cell data with clustering done in h5ad format.

  	output_name
  		Output file name. All outputs will use it as the prefix.

* Options:

	-\\-subset-selection <subset-selection>...
		Specify which cells will be included in the subcluster analysis. Each <subset_selection> string takes the format of 'attr:value,...,value', which means select cells with attr in the values. If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings.

	\-p <number>, -\\-threads <number>
		Number of threads. [default: 1]

	-\\-correct-batch-effect
		Correct for batch effects.

	-\\-batch-group-by
		Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+pegasus..+attrn', or 'attr=value11,pegasus..,value1n_1;value21,pegasus..,value2n_2;pegasus..;valuem1,pegasus..,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',pegasus..,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,pegasus..,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.

	-\\-output-loom
		Output loom-formatted file.

	-\\-select-hvf-flavor <flavor>
		Highly variable feature selection method. <flavor> can be 'pegasus' or 'Seurat'. [default: pegasus]

	-\\-select-hvf-ngenes <nfeatures>
		Select top <nfeatures> highly variable features. If <flavor> is 'Seurat' and <nfeatures> is 'None', select HVGs with z-score cutoff at 0.5 [default: 2000]

	-\\-no-select-hvf
		Do not select highly variable features.

	-\\-plot-hvf
		Plot highly variable feature selection.

	-\\-random-state <seed>
		Random number generator seed. [default: 0]

	-\\-temp-folder <temp_folder>
		Joblib temporary folder for memmapping numpy arrays.
  
	-\\-nPC <number>
		Number of principal components. [default: 50]

	-\\-knn-K <number>
		Number of nearest neighbors for building kNN graph. [default: 100]

	-\\-knn-full-speed
		For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

	-\\-kBET
		Calculate kBET.

	-\\-kBET-batch <batch>
		kBET batch keyword.

	-\\-kBET-alpha <alpha>
		kBET rejection alpha. [default: 0.05]

	-\\-kBET-K <K> 
		kBET K. [default: 25]

	-\\-diffmap
		Calculate diffusion maps.

	-\\-diffmap-ndc <number>
		Number of diffusion components. [default: 50]

	-\\-diffmap-alpha <alpha>
		Power parameter for diffusion-based pseudotime. [default: 0.5]

	-\\-diffmap-solver <solver>
		Solver for eigen decomposition, either 'randomized' or 'eigsh'. [default: randomized]

	-\\-diffmap-to-3d
		If map diffusion map into 3D space using PCA.

	-\\-calculate-pseudotime <roots>
		Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  	-\\-louvain
  		Run louvain clustering algorithm.

	-\\-louvain-resolution <resolution>
		Resolution parameter for the louvain clustering algorithm. [default: 1.3]

	-\\-louvain-class-label <label>
		Louvain cluster label name in AnnData. [default: louvain_labels]

	-\\-leiden
		Run leiden clustering algorithm.

	-\\-leiden-resolution <resolution>
		Resolution parameter for the leiden clustering algorithm. [default: 1.3]

	-\\-leiden-niter <niter>
		Number of iterations of running the Leiden algorithm. If <niter> is negative, run Leiden iteratively until no improvement. [default: -1]

	-\\-leiden-class-label <label>
		Leiden cluster label name in AnnData. [default: leiden_labels]

	-\\-spectral-louvain
		Run spectral-louvain clustering algorithm.

	-\\-spectral-louvain-basis <basis>
		Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. [default: diffmap]

	-\\-spectral-louvain-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-spectral-louvain-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-spectral-louvain-resolution <resolution>.
		Resolution parameter for louvain. [default: 1.3]

	-\\-spectral-louvain-class-label <label>
		Spectral-louvain label name in AnnData. [default: spectral_louvain_labels]

	-\\-spectral-leiden
		Run spectral-leiden clustering algorithm.

	-\\-spectral-leiden-basis <basis>
		Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. [default: diffmap]

	-\\-spectral-leiden-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-spectral-leiden-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-spectral-leiden-resolution <resolution>
		Resolution parameter for leiden. [default: 1.3]

	-\\-spectral-leiden-class-label <label>
		Spectral-leiden label name in AnnData. [default: spectral_leiden_labels]

	-\\-tsne
		Run multi-core t-SNE for visualization.

	-\\-run-fitsne
  		Run FIt-SNE for visualization.

	-\\-tsne-perplexity <perplexity>
		t-SNE's perplexity parameter. [default: 30]

  	-\\-umap
  		Run umap for visualization.

	-\\-umap-K <K>
		K neighbors for umap. [default: 15]

	-\\-umap-min-dist <number>
		Umap parameter. [default: 0.5]

	-\\-umap-spread <spread>
		Umap parameter. [default: 1.0]

	-\\-fle
		Run force-directed layout embedding.

	-\\-fle-K <K>
		K neighbors for building graph for FLE. [default: 50]

	-\\-fle-target-change-per-node <change>
		Target change per node to stop forceAtlas2. [default: 2.0]

	-\\-fle-target-steps <steps>
		Maximum number of iterations before stopping the forceAtlas2 algoritm. [default: 5000]

	-\\-fle-memory <memory>
		Memory size in GB for the Java FA2 component. [default: 8]

	-\\-net-down-sample-fraction <frac>
		Down sampling fraction for net-related visualization. [default: 0.1]

	-\\-net-down-sample-K <K>
		Use <K> neighbors to estimate local density for each data point for down sampling. [default: 25]

	-\\-net-down-sample-alpha <alpha>
		 Weighted down sample, proportional to radius^alpha. [default: 1.0]

	-\\-net-regressor-L2-penalty <value>
		L2 penalty parameter for the deep net regressor. [default: 0.1]

	-\\-net-tsne
		Run net tSNE for visualization.

	-\\-net-tsne-polish-learning-frac <frac>
		After running the deep regressor to predict new coordinates, use <frac> * nsample as the learning rate to use to polish the coordinates. [default: 0.33]

	-\\-net-tsne-polish-niter <niter>
		Number of iterations for polishing tSNE run. [default: 150]

	-\\-net-tsne-out-basis <basis>
		Output basis for net-tSNE. [default: net_tsne]

	-\\-net-umap
		Run net umap for visualization.

	-\\-net-umap-polish-learning-rate <rate>
		After running the deep regressor to predict new coordinate, what is the learning rate to use to polish the coordinates for UMAP. [default: 1.0]

	-\\-net-umap-polish-nepochs <nepochs>
		Number of iterations for polishing UMAP run. [default: 40]

	-\\-net-umap-out-basis <basis>
		Output basis for net-UMAP. [default: net_umap]

	-\\-net-fle
		Run net FLE.

	-\\-net-fle-polish-target-steps <steps>
		After running the deep regressor to predict new coordinate, what is the number of force atlas 2 iterations. [default: 1500]

	-\\-net-fle-out-basis <basis>
		Output basis for net-FLE. [default: net_fle]

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.h5ad
		Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, t-SNE and diffusion map coordinates are stored in the 'obsm' field.

	output_name.loom
		Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.

* Examples::

	pegasus subcluster -p 20 --correct-batch-effect --subset-selection louvain_labels:3,6 --subset-selection Condition:CB_nonmix --tsne --louvain manton_bm.h5ad manton_bm_subset


---------------------------------


``pegasus scp_output``
^^^^^^^^^^^^^^^^^^^^^^^

If we want to visualize analysis results on single cell portal (SCP), we can generate required files for SCP using this subcommand.

Type::

	pegasus scp_output -h

to see the usage information::

	Usage:
		pegasus scp_output <input_h5ad_file> <output_name>
		pegasus scp_output -h

* Arguments:

	input_h5ad_file
		Analyzed single cell data in h5ad format.

	output_name
		Name prefix for all outputted files.

* Options:

	-\\-dense
		Output dense expression matrix instead.

	-\\-round-to <ndigit>
		Round expression to <ndigit> after the decimal point. [default: 2]

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.scp.metadata.txt, output_name.scp.barcodes.tsv, output_name.scp.genes.tsv, output_name.scp.matrix.mtx, output_name.scp.*.coords.txt
		Files that single cell portal needs.

* Examples::

	pegasus scp_output example.h5ad example


---------------------------------


``pegasus parquet``
^^^^^^^^^^^^^^^^^^^^^^^

Generate a PARQUET file for web-based visualization.

Type::

	pegasus parquet -h

to see the usage information::

	Usage:
		pegasus parquet [options] <input_h5ad_file> <output_name>
		pegasus parquet -h

* Arguments:

	input_h5ad_file
		Analyzed single cell data in h5ad format.

	output_name
		Name prefix for the parquet file.

* Options:

	\-p <number>, -\\-threads <number>
		Number of threads used to generate the PARQUET file. [default: 1]

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.parquet
		Generated PARQUET file that contains metadata and expression levels for every gene.

* Examples::

	pegasus parquet example.h5ad example.parquet


---------------------------------


``pegasus merge_rna_adt``
^^^^^^^^^^^^^^^^^^^^^^^^^

If we have CITE-Seq data, we can merge RNA count matrix and ADT (antibody tag) count matrix into one file using this subcommand.

Type::

	pegasus merge_rna_adt -h

to see the usage information::

	Usage:
		pegasus merge_rna_adt <input_raw_gene_bc_matrices_h5.h5sc> <input_adt_csv_file> <output_name>
		pegasus merge_rna_adt -h

* Arguments:

	input_raw_gene_bc_matrices_h5.h5sc
		Input raw RNA expression matrix in pegasus hdf5 format.

	input_adt_csv_file
		Input ADT (antibody tag) count matrix in CSV format.

	output_name
		Merged output name.

* Options:

	-\\-antibody-control-csv <antibody_control_csv_file>
		A CSV file containing the IgG control information for each antibody.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.h5sc
		Output file in pegasus hdf5 format. This file contains two groups --- one for RNAs, and the other for ADTs.

* Examples::

	pegasus merge_rna_adt example_raw_h5.h5sc example_adt.csv example_merged_raw
	pegasus merge_rna_adt --antibody-control-csv antibody_control.csv example_raw_h5.h5sc example_adt.csv example_merged_raw


---------------------------------


``pegasus check_indexes``
^^^^^^^^^^^^^^^^^^^^^^^^^

If we run CITE-Seq or any kind of hashing, we need to make sure that the library indexes of CITE-Seq/hashing do not collide with 10x's RNA indexes. This command can help us to determine which 10x index sets we should use.

Type::

	pegasus check_indexes -h

to see the usage information::

	Usage:
		pegasus check_indexes [--num-mismatch <mismatch> --num-report <report>] <index_file>
		pegasus check_indexes -h

* Arguments:

	index_file
		Index file containing CITE-Seq/hashing index sequences. One sequence per line.

* Options:

	-\\-num-mismatch <mismatch>
		Number of mismatch allowed for each index sequence. [default: 1]

  	-\\-num-report <report>
  		Number of valid 10x indexes to report. Default is to report all valid indexes. [default: 9999]
  
  	\-h, -\\-help
  		Print out help information.

* Outputs:

	Up to <report> number of valid 10x indexes will be printed out to standard output.

* Examples::

	pegasus check_indexes --num-report 8 index_file.txt

