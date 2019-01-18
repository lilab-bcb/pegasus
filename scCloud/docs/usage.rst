Use ``scCloud`` as a command line tool
---------------------------------------

``scCloud`` can be used as a command line tool. Type::

	scCloud -h

to see the help information::

	Usage:
		scCloud <command> [<args>...]
		scCloud -h | --help
		scCloud -v | --version

``scCloud`` has 13 sub-commands in 8 groups.

* Preprocessing:

	aggregate_matrix
		Aggregate cellranger-outputted channel-specific count matrices into a single count matrix. It also enables users to import metadata into the count matrix.

* Demultiplexing:

	demuxEM
		Demultiplex cells/nuclei based on DNA barcodes for cell-hashing and nuclei-hashing data.

* Analyzing:
	
	cluster
		Perform first-pass analysis using the count matrix generated from 'aggregate_matrix'. This subcommand could perform low quality cell filtration, batch correction, variable gene selection, dimension reduction, diffusion map calculation, graph-based clustering, tSNE visualization. The final results will be written into h5ad-formatted file, which Seurat could load.
  		
    de_analysis
    	Detect markers for each cluster by performing differential expression analysis per cluster (within cluster vs. outside cluster). DE tests include Welch's t-test, Fisher's exact test, Mann-Whitney U test. It can also calculate AUROC values for each gene.
    
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

	scCloud aggregate_matrix --attributes Source,Platform,Donor example.csv example
	scCloud cluster -p 20 --correct-batch-effect --batch-group-by Source -run-louvain --run-tsne example_10x.h5 example
	scCloud de_analysis --labels louvain_labels -p 20 --fisher example.h5ad example_de.xlsx
	scCloud annotate_cluster example.h5ad example.anno.txt
	scCloud plot composition --cluster-labels louvain_labels --attribute Donor --style normalized --not-stacked example.h5ad example.composition.pdf
	scCloud plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad example.scatter.pdf
	scCloud iplot --attribute louvain_labels diffmap_pca example.h5ad example.diffmap.html

The above analysis will give you tSNE, louvain cluster labels and diffusion maps in ``example.h5ad``. You can investigate donor-specific effects by looking at ``example.composition.pdf``. ``example.scatter.pdf`` plotted tSNE colored by louvain_labels and Donor info side-by-side. You can explore the diffusion map in 3D by looking at ``example.diffmap.html``. This html maps all diffusion components into 3D using PCA.

If you want to perform subcluster analysis by combining cluster 1 and 3, run the following command::

	scCloud subcluster -p 20 --correct-batch-effect example.h5ad 1,3 example_sub


---------------------------------


``scCloud aggregate_matrix``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step for single cell analysis is to generate one count matrix from cellranger's channel-specific count matrices. ``scCloud aggregate_matrix`` allows aggregating arbitrary matrices with the help of a *CSV* file.

Type::

	scCloud aggregate_matrix -h

to see the usage information::

	Usage:
		scCloud aggregate_matrix <csv_file> <output_name> [--restriction <restriction>... --attributes <attributes> --google-cloud --select-only-singlets --minimum-number-of-genes <ngene> --dropseq-genome <genome>]
		scCloud aggregate_matrix -h

* Arguments:

	csv_file
		Input csv-formatted file containing information of each scRNA-Seq run. Each row must contain at least 2 columns --- Sample, sample name and Location, location of the channel-specific count matrix in either 10x format (e.g. /sample/filtered_gene_bc_matrices_h5.h5) or dropseq format (e.g. /sample/sample.umi.dge.txt.gz). See below for an example csv::

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

	-\\-google-cloud
		If files are stored in google cloud. Assuming google cloud sdk is installed.

	-\\-select-only-singlets
		If we have demultiplexed data, turning on this option will make scCloud only include barcodes that are predicted as singlets.

	-\\-minimum-number-of-genes <ngene>
		Only keep barcodes with at least <ngene> expressed genes.

	-\\-dropseq-genome <genome>
		If inputs are dropseq data, this option needs to turn on and provides the reference genome name.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name_10x.h5
		A 10x-formatted HDF5 file containing the count matrices and associated attributes.

* Examples::

	scCloud aggregate_matrix --restriction Source:pbmc --restriction Donor:1 --attributes Source,Platform,Donor example.csv example


---------------------------------


``scCloud demuxEM``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have data generated by cell-hashing or nuclei-hashing, you can use ``scCloud demuxEM`` to demultiplex your data. 

Type::

	scCloud demuxEM -h

to see the usage information::

	Usage:
		scCloud demuxEM --hash-type <type> [options] <input_adt_csv_file> <input_raw_gene_bc_matrices_h5.h5> <output_name>
		scCloud demuxEM -h

* Arguments:

	input_adt_csv_file
		Input ADT (antibody tag) count matrix in CSV format.

	input_raw_gene_bc_matrices_h5.h5
		Input raw RNA expression matrix in 10x hdf5 format.

	output_name
		Output name. All outputs will use it as the prefix.

* Options:

	-\\-hash-type <type>
		The hash type of the data. <type> can be 'cell-hashing' for cell-hashing and 'nuclei-hashing' for nuclei-hashing.

  	\-p <number>, -\\-threads <number>
		Number of threads. [default: 1]

	-\\-genome <genome>
		Reference genome name. If not provided, we will infer it from the expression matrix file.

	-\\-min-num-genes <number>
		We only demultiplex cells/nuclei with at least <number> expressed genes. [default: 100]

	-\\-min-signal-hashtag <count>
		Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown. [default: 10.0]

	-\\-prior-on-samples <prior>
		The sparse prior put on samples.

	-\\-random-state <seed>
		The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]

	-\\-generate-diagnostic-plots
		Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc.

	-\\-generate-gender-plot <genes>
		Generate violin plots using gender-specific genes (e.g. Xist). <gene> is a comma-separated list of gene names. 
	
	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name_demux_10x.h5
		RNA expression matrix with demultiplexed sample identities in 10x's hdf5 format.

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

	scCloud demuxEM -p 8 --hash-type cell-hashing --generate-diagnostic-plots example_adt.csv example_raw_gene_bc_matrices_h5.h5 example_output


---------------------------------


``scCloud cluster``
^^^^^^^^^^^^^^^^^^^

Once we collected the count matrix ``example_10x.h5``, we can perform single cell analysis using ``scCloud cluster``.

Type::

	scCloud cluster -h

to see the usage information::

	Usage:
		scCloud cluster [options] <input_file> <output_name>
		scCloud cluster -h

* Arguments:

	input_file
		Input file in 10x format. If first-pass analysis has been performed, but you want to run some additional analysis, you could also pass a h5ad-formatted file.

	output_name      
		Output file name. All outputs will use it as the prefix.

* Options:

	\-p <number>, -\\-threads <number>
		Number of threads. [default: 1]

	-\\-processed
		Input file is processed and thus no PCA & diffmap will be run.

	-\\-genome <genome>
		A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
  
	-\\-cite-seq
		Data are CITE-Seq data. scCloud will perform analyses on RNA count matrix first. Then it will attach the ADT matrix to the RNA matrix with all antibody names changing to 'AD-' + antibody_name. Lastly, it will embed the antibody expression using t-SNE (the basis used for plotting is 'citeseq_tsne').

  	-\\-output-filtration-results
		Output filtration results as a spreadsheet.

	-\\-plot-filtration-results
		Plot filtration results as PDF files.

	-\\-plot-filtration-figsize <figsize>
		Figure size for filtration plots. <figsize> is a comma-separated list of two numbers, the width and height of the figure (e.g. 6,4).

	-\\-make-output-seurat-compatible
		Make output h5ad file seurat compatible. Caution: this will significantly increase the output size. Do not turn this option on for large data sets.

	-\\-output-loom
		Output loom-formatted file.

	-\\-correct-batch-effect
		Correct for batch effects.

	-\\-batch-group-by <expression>
		Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+...+attrn', or 'attr=value11,...,value1n_1;value21,...,value2n_2;...;valuem1,...,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',...,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,...,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.

  	-\\-min-genes <number>
		Only keep cells with at least <number> of genes. [default: 500]

	-\\-max-genes <number>
		Only keep cells with less than <number> of genes. [default: 6000]

	-\\-min-umis <number>
		Only keep cells with at least <number> of UMIs. [default: 100]

	-\\-max-umis <number>
		Only keep cells with less than <number> of UMIs. [default: 600000]

	-\\-mito-prefix <prefix>
		Prefix for mitochondrial genes. [default: MT-]

	-\\-percent-mito <ratio>
		Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]

	-\\-gene-percent-cells <ratio>
		Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]

	-\\-min-genes-on-raw <number>
		If input are raw 10x matrix, which include all barcodes, perform a pre-filtration step to keep the data size small. In the pre-filtration step, only keep cells with at least <number> of genes. [default: 100]

	-\\-counts-per-cell-after <number>
		Total counts per cell after normalization. [default: 1e5]

	-\\-random-state <seed>
		Random number generator seed. [default: 0]

	-\\-run-uncentered-pca
		Run uncentered PCA.

	-\\-no-variable-gene-selection
		Do not select variable genes.

	-\\-no-submat-to-dense
		Do not convert variable-gene-selected submatrix to a dense matrix.
  
	-\\-nPC <number>
		Number of PCs. [default: 50]

	-\\-nDC <number>
		Number of diffusion components. [default: 50]

	-\\-diffmap-alpha <alpha>
		Power parameter for diffusion-based pseudotime. [default: 0.5]

	-\\-diffmap-K <K>
		Number of neighbors used for constructing affinity matrix. [default: 100]

	-\\-diffmap-full-speed
		For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

	-\\-calculate-pseudotime <roots>
		Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  	-\\-run-louvain
  		Run louvain clustering algorithm.

	-\\-louvain-resolution <resolution>
		Resolution parameter for the louvain clustering algorithm. [default: 1.3]

	-\\-louvain-affinity <affinity>
		Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

	-\\-run-approximated-louvain
		Run approximated louvain clustering algorithm.

	-\\-approx-louvain-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-approx-louvain-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-approx-louvain-resolution <resolution>.
		Resolution parameter for louvain. [default: 1.3]

	-\\-run-tsne
		Run multi-core t-SNE for visualization.

	-\\-tsne-perplexity <perplexity>
		t-SNE's perplexity parameter. [default: 30]

  	-\\-run-fitsne
  		Run FIt-SNE for visualization.

  	-\\-run-umap
  		Run umap for visualization.

	-\\-umap-on-diffmap
		Run umap on diffusion components.

	-\\-umap-K <K>
		K neighbors for umap. [default: 15]

	-\\-umap-min-dist <number>
		Umap parameter. [default: 0.1]

	-\\-umap-spread <spread>
		Umap parameter. [default: 1.0]

	-\\-run-fle
		Run force-directed layout embedding.

	-\\-fle-K <K>
		K neighbors for building graph for FLE. [default: 50]

	-\\-fle-n-steps <nstep>
		Number of iterations for FLE. [default: 10000]

	-\\-fle-affinity <affinity>
		Affinity matrix to be used. Could be 'W_diffmap', or 'W_diffmap_norm'. [default: W_diffmap]

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.h5ad
		Output file in h5ad format. To load this file in python, use ``import scCloud; data = scCloud.tools.read_input('output_name.h5ad', mode = 'a')``. The log-normalized expression matrix is stored in ``data.X`` as a CSR-format sparse matrix. The ``obs`` field contains cell related attributes, including clustering results. For example, ``data.obs_names`` records cell barcodes; ``data.obs['Channel']`` records the channel each cell comes from; ``data.obs['n_genes']``, ``data.obs['n_counts']``, and ``data.obs['percent_mito']`` record the number of expressed genes, total UMI count, and mitochondrial rate for each cell respectively; ``data.obs['louvain_labels']`` and ``data.obs['approx_louvain_labels']`` record each cell's cluster labels using different clustring algorithms; ``data.obs['pseudo_time']`` records the inferred pseudotime for each cell. The ``var`` field contains gene related attributes. For example, ``data.var_names`` records gene symbols, ``data.var['gene_ids']`` records Ensembl gene IDs, and ``data.var['selected']`` records selected variable genes. The ``obsm`` field records embedding coordiates. For example, ``data.obsm['X_pca']`` records PCA coordinates, ``data.obsm['X_tsne']`` records tSNE coordinates, ``data.obsm['X_umap']`` records UMAP coordinates, ``data.obsm['X_diffmap']`` records diffusion map coordinates, ``data.obsm['X_diffmap_pca']`` records the first 3 PCs by projecting the diffusion components using PCA, and ``data.obsm['X_fle']`` records the force-directed layout coordinates from the diffusion components. The ``uns`` field stores other related information, such as reference genome (``data.uns['genome']``). If '--make-output-seurat-compatible' is on, this file can be loaded into R and converted into a Seurat object.

	output_name.filt.xlsx
		Optional output. Only exists if '--output-filtration-results' is set. This file has two sheets --- Cell filtration stats and Gene filtration stats. The first sheet records cell filtering results and it has 10 columns: Channel, channel name; kept, number of cells kept; median_n_genes, median number of expressed genes in kept cells; median_n_umis, median number of UMIs in kept cells; median_percent_mito, median mitochondrial rate as UMIs between mitochondrial genes and all genes in kept cells; filt, number of cells filtered out; total, total number of cells before filtration, if the input contain all barcodes, this number is the cells left after '--min-genes-on-raw' filtration; median_n_genes_before, median expressed genes per cell before filtration; median_n_umis_before, median UMIs per cell before filtration; median_percent_mito_before, median mitochondrial rate per cell before filtration. The channels are sorted in ascending order with respect to the number of kept cells per channel. The second sheet records genes that failed to pass the filtering. This sheet has 3 columns: gene, gene name; n_cells, number of cells this gene is expressed; percent_cells, the fraction of cells this gene is expressed. Genes are ranked in ascending order according to number of cells the gene is expressed. Note that only genes not expressed in any cell are removed from the data. Other filtered genes are marked as non-robust and not used for TPM-like normalization.

	output_name.filt.gene.pdf
		Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting gene count distributions before and after filtration per channel.
	
	output_name.filt.UMI.pdf
		Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting UMI count distributions before and after filtration per channel.
	
	output_name.filt.mito.pdf
		Optional output. Only exists if '--plot-filtration-results' is set. This file contains violin plots contrasting mitochondrial rate distributions before and after filtration per channel.

	output_name.loom
		Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.

* Examples::

	scCloud cluster -p 20 --correct-batch-effect --run-louvain --run-tsne example_10x.h5 example


---------------------------------


``scCloud de_analysis``
^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the clusters, we can detect markers using ``scCloud de_analysis``.

Type::

	scCloud de_analysis -h

to see the usage information::

	Usage:
		scCloud de_analysis [--labels <attr> -p <threads> --alpha <alpha> --fisher --mwu --roc] <input_h5ad_file> <output_spreadsheet>
		scCloud de_analysis -h

* Arguments:

	input_h5ad_file
		Single cell data with clustering calculated. DE results would be written back.
	
	output_spreadsheet
		Output spreadsheet with DE results.

* Options:

	-\\-labels <attr>
		<attr> used as cluster labels. [default: louvain_labels]

	-\\-alpha <alpha>
		Control false discovery rate at <alpha>. [default: 0.05]

	-\\-fisher
		Calculate Fisher's exact test.

	-\\-mwu
		Calculate Mann-Whitney U test.

	-\\-roc
		Calculate area under cuver in ROC curve.

	\-p <threads>
		Use <threads> threads. [default: 1]

	\-h, -\\-help
		Print out help information.

* Outputs:

	input_h5ad_file
		DE results would be written back to the 'var' fields.

	output_spreadsheet
		An excel spreadsheet containing DE results. Each cluster has two tabs in the spreadsheet. One is for up-regulated genes and the other is for down-regulated genes.

* Examples::

	scCloud de_analysis --labels louvain_labels -p 20 --fisher --mwu --roc example.h5ad example_de.xlsx


---------------------------------


``scCloud annotate_cluster``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the DE results, we could optionally identify putative cell types for each cluster using ``scCloud annotate_cluster``. Currently, this subcommand works for human/mouse immune/brain cells. This command has two forms: the first form generates putative annotations and the second form write annotations into the h5ad object.

Type::

	scCloud annotate_cluster -h

to see the usage information::

	Usage:
		scCloud annotate_cluster [--json-file <file> --minimum-report-score <score> --do-not-use-non-de-genes] <input_h5ad_file> <output_file>
		scCloud annotate_cluster --annotation <annotation_string> <input_h5ad_file>
		scCloud annotate_cluster -h

* Arguments:

	input_h5ad_file
		Single cell data with DE analysis done by ``scCloud de_analysis``.

	output_file
		Output annotation file.

* Options:

	-\\-json-file <file>
		JSON file for markers. Could also be ``human_immune``/``mouse_immune``/``mouse_brain``/``human_brain``, which triggers scCloud to markers included in the package. [default: human_immune]

	-\\-minimum-report-score <score>
		Minimum cell type score to report a potential cell type. [default: 0.5]

	-\\-do-not-use-non-de-genes
		Do not count non DE genes as down-regulated.

	-\\-annotation <annotation_string>
		Write cell type annotations in <annotation_string> into <input_h5ad_file>. <annotation_string> has this format: 'anno_attr:anno_1;anno_2;...;anno_n'. 'anno_attr' is the annotation attribute in the h5ad object and anno_i is the annotation for cluster i.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_file
		This is a text file. For each cluster, all its putative cell types are listed in descending order of the cell type score. For each putative cell type, all markers support this cell type are listed. If one putative cell type has cell subtypes, all subtypes will be listed under this cell type.

* Examples::

	scCloud annotate_cluster example.h5ad example.anno.txt
	scCloud annotate_cluster --annotation "anno:T cells;B cells;NK cells;Monocytes" example.h5ad


---------------------------------



``scCloud plot``
^^^^^^^^^^^^^^^^^

We can make a variety of figures using ``scCloud plot``.

Type::

	scCloud plot -h

to see the usage information::

	Usage:
  		scCloud plot [options] [--restriction <restriction>...] <plot_type> <input_h5ad_file> <output_file>
		scCloud plot -h

* Arguments:

	plot_type
		Only 2D plots, chosen from 'composition', 'scatter', 'scatter_groups', 'scatter_genes', 'scatter_gene_groups', and 'heatmap'.

	input_h5ad_file
		Single cell data with clustering done by Scanpy in h5ad file format.

  	output_file
  		Output image file.

* Options:

	-\\-dpi <dpi>
		DPI value for the figure. [default: 500]

	-\\-cluster-labels <attr>
		Use <attr> as cluster labels. This option is used in 'composition', 'scatter_groups', and 'heatmap'.

  	-\\-attribute <attr>
  		Plot <attr> against cluster labels. This option is only used in 'composition'.

	-\\-basis <basis>
		Basis for 2D plotting, chosen from 'tsne', 'fitsne', 'umap', 'pca', 'rpca', 'fle', or 'diffmap_pca'. If CITE-Seq data is used, basis can also be 'citeseq_tsne'. This option is used in 'scatter', 'scatter_groups', 'scatter_genes', and 'scatter_gene_groups'. [default: tsne]

	-\\-attributes <attrs>
		<attrs> is a comma-separated list of attributes to color the basis. This option is only used in 'scatter'.

	-\\-restriction <restriction>...
		Set restriction if you only want to plot a subset of data. Multiple <restriction> strings are allowed. Each <restriction> takes the format of 'attr:value,value'. This option is used in 'composition' and 'scatter'.
	
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
		Number of rows in the figure. If not set, scCloud will figure it out automatically.

	-\\-ncols <ncols>
		Number of columns in the figure. If not set, scCloud will figure it out automatically.

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

	\-h, -\\-help
		Print out help information.

Examples::

	scCloud plot composition --cluster-labels louvain_labels --attribute Donor --style normalized --not-stacked example.h5ad example.composition.pdf
	scCloud plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad example.scatter.pdf
	scCloud plot scatter_groups --cluster-labels louvain_labels --group Donor example.h5ad example.scatter_groups.pdf
	scCloud plot scatter_genes --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP example.h5ad example.genes.pdf
	scCloud plot scatter_gene_groups --gene CD8A --group Donor example.h5ad example.gene_groups.pdf
	scCloud plot heatmap --cluster-labels louvain_labels --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP --heatmap-title 'markers' example.h5ad example.heatmap.pdf


---------------------------------


``scCloud iplot``
^^^^^^^^^^^^^^^^^^

We can also make interactive plots in html format using ``scCloud iplot``. These interactive plots are very helpful if you want to explore the diffusion maps.

Type::

	scCloud iplot -h

to see the usage information::

	Usage:
		scCloud iplot --attribute <attr> [options] <basis> <input_h5ad_file> <output_html_file>
		scCloud iplot -h

* Arguments:

	basis
		Basis can be either 'tsne', 'fitsne', 'umap', 'diffmap', 'pca', 'rpca' or 'diffmap_pca'.
	
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

	scCloud iplot --attribute louvain_labels tsne example.h5ad example.tsne.html
	scCloud iplot --attribute louvain_labels diffmap_pca example.h5ad example.diffmap.html


---------------------------------


``scCloud view``
^^^^^^^^^^^^^^^^^

We may want to further perform sub-cluster analysis on a subset of cells. This sub-command helps us to define the subset.

Type::

	scCloud view -h

to see the usage information::

	Usage:
		scCloud view [--show-attributes --show-gene-attributes --show-values-for-attributes <attributes>] <input_h5ad_file>
		scCloud view -h

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

	scCloud view --show-attributes example.h5ad
	scCloud view --show-gene-attributes example.h5ad
	scCloud view --show-values-for-attributes louvain_labels,Donor example.h5ad


---------------------------------


``scCloud subcluster``
^^^^^^^^^^^^^^^^^^^^^^^

If there is a subset of cells that we want to further cluster, we can run ``scCloud subcluster``. This sub-command will outputs a new h5ad file that you can run ``de_analysis``, ``plot`` and ``iplot`` on.

Type::

	scCloud subcluster -h

to see the usage information::

	Usage:
		scCloud subcluster [options] --subset-selection <subset-selection>... <input_file> <output_name>
		scCloud subcluster -h

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

	-\\-output-loom
		Output loom-formatted file.

	-\\-random-state <seed>
		Random number generator seed. [default: 0]

	-\\-run-uncentered-pca
		Run uncentered PCA.

	-\\-no-variable-gene-selection
		Do not select variable genes.

	-\\-no-submat-to-dense
		Do not convert variable-gene-selected submatrix to a dense matrix.
  
	-\\-nPC <number>
		Number of PCs. [default: 50]

	-\\-nDC <number>
		Number of diffusion components. [default: 50]

	-\\-diffmap-alpha <alpha>
		Power parameter for diffusion-based pseudotime. [default: 0.5]

	-\\-diffmap-K <K>
		Number of neighbors used for constructing affinity matrix. [default: 100]

	-\\-diffmap-full-speed
		For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.

	-\\-calculate-pseudotime <roots>
		Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.

  	-\\-run-louvain
  		Run louvain clustering algorithm.

	-\\-louvain-resolution <resolution>
		Resolution parameter for the louvain clustering algorithm. [default: 1.3]

	-\\-louvain-affinity <affinity>
		Affinity matrix to be used. Could be 'W_norm', 'W_diffmap', or 'W_diffmap_norm'. [default: W_norm]

	-\\-run-approximated-louvain
		Run approximated louvain clustering algorithm.

	-\\-approx-louvain-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-approx-louvain-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-approx-louvain-resolution <resolution>.
		Resolution parameter for louvain. [default: 1.3]

	-\\-run-tsne
		Run multi-core t-SNE for visualization.

	-\\-tsne-perplexity <perplexity>
		t-SNE's perplexity parameter. [default: 30]

  	-\\-run-fitsne
  		Run FIt-SNE for visualization.

  	-\\-run-umap
  		Run umap for visualization.

	-\\-umap-on-diffmap
		Run umap on diffusion components.

	-\\-umap-K <K>
		K neighbors for umap. [default: 15]

	-\\-umap-min-dist <number>
		Umap parameter. [default: 0.1]

	-\\-umap-spread <spread>
		Umap parameter. [default: 1.0]

	-\\-run-fle
		Run force-directed layout embedding.

	-\\-fle-K <K>
		K neighbors for building graph for FLE. [default: 50]

	-\\-fle-n-steps <nstep>
		Number of iterations for FLE. [default: 10000]

	-\\-fle-affinity <affinity>
		Affinity matrix to be used. Could be 'W_diffmap', or 'W_diffmap_norm'. [default: W_diffmap]

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.h5ad
		Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, t-SNE and diffusion map coordinates are stored in the 'obsm' field.

	output_name.loom
		Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.

* Examples::

	scCloud subcluster --subset_selection louvain_labels:1,3  --subset_selection Donor:1 -p 20 --correct-batch-effect example.h5ad example_sub


---------------------------------


``scCloud scp_output``
^^^^^^^^^^^^^^^^^^^^^^^

If we want to visualize analysis results on single cell portal (SCP), we can generate required files for SCP using this subcommand.

Type::

	scCloud scp_output -h

to see the usage information::

	Usage:
		scCloud scp_output <input_h5ad_file> <output_name>
		scCloud scp_output -h

* Arguments:

	input_h5ad_file
		Analyzed single cell data in h5ad format.

	output_name
		Name prefix for all outputted files.

* Options:

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name.scp.metadata.txt, output_name.scp.barcodes.tsv, output_name.scp.genes.tsv, output_name.scp.matrix.mtx, output_name.scp.*.coords.txt
		Files that single cell portal needs.

* Examples::

	scCloud scp_output example.h5ad example


---------------------------------


``scCloud parquet``
^^^^^^^^^^^^^^^^^^^^^^^

Generate a PARQUET file for web-based visualization.

Type::

	scCloud parquet -h

to see the usage information::

	Usage:
		scCloud parquet [options] <input_h5ad_file> <output_name>
		scCloud parquet -h

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

	scCloud parquet manton_bm.h5ad manton_bm.parquet


---------------------------------


``scCloud merge_rna_adt``
^^^^^^^^^^^^^^^^^^^^^^^^^

If we have CITE-Seq data, we can merge RNA count matrix and ADT (antibody tag) count matrix into one file using this subcommand.

Type::

	scCloud merge_rna_adt -h

to see the usage information::

	Usage:
		scCloud merge_rna_adt <input_raw_gene_bc_matrices_h5.h5> <input_adt_csv_file> <antibody_control_csv> <output_10x.h5>
		scCloud merge_rna_adt -h

* Arguments:

	input_raw_gene_bc_matrices_h5.h5
		Input raw RNA expression matrix in 10x hdf5 format.

	input_adt_csv_file
		Input ADT (antibody tag) count matrix in CSV format.

	antibody_control_csv
		A CSV file containing the IgG control information for each antibody.

	output_10x.h5
		Merged output file in 10x hdf5 format.

* Options:

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_10x.h5
		Output file in 10x hdf5 format. This file contains two groups --- one is for RNAs and the other is for ADTs.

* Examples::

	scCloud merge_rna_adt example_raw_h5.h5 example_adt.csv antibody_control.csv example_merged_raw_10x.h5


---------------------------------


``scCloud check_indexes``
^^^^^^^^^^^^^^^^^^^^^^^^^

If we run CITE-Seq or any kind of hashing, we need to make sure that the library indexes of CITE-Seq/hashing do not collide with 10x's RNA indexes. This command can help us to determine which 10x index sets we should use.

Type::

	scCloud check_indexes -h

to see the usage information::

	Usage:
		scCloud check_indexes [--num-mismatch <mismatch> --num-report <report>] <index_file>
		scCloud check_indexes -h

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

	scCloud check_indexes --num-report 8 index_file.txt

