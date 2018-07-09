Use ``scrtools`` as a command line tool
---------------------------------------

``scrtools`` can be used as a command line tool. Type::

	scrtools -h

to see the help information::

	Usage:
		scrtools <command> [<args>...]
		scrtools -h | --help
		scrtools -v | --version

``scrtools`` has 8 sub-commands in 4 groups.

* Preprocessing:

	aggregate_matrix
		Aggregate cellranger-outputted channel-specific count matrices into a single count matrix. It also enables users to import metadata into the count matrix.

* Analyzing:
	
	cluster
		Perform first-pass analysis using the count matrix generated from 'aggregate_matrix'. This subcommand could perform low quality cell filtration, batch correction, variable gene selection, dimension reduction, diffusion map calculation, graph-based clustering, tSNE visualization. The final results will be written into h5ad-formatted file, which Seurat could load.
  		
    de_analysis
    	Detect markers for each cluster by performing differential expression analysis per cluster (within cluster vs. outside cluster). DE tests include Welch's t-test, Fisher's exact test, Mann-Whitney U test. It can also calculate AUROC values for each gene.
    
    annotate_cluster
    	This subcommand is used to automatically annotate cell types for each cluster based on existing markers. Currently, it only works for human and mouse immune cells.

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

	scrtools aggregate_matrix --genome GRCh38 --attributes Source,Platform,Donor example.csv example
	scrtools cluster -p 20 --correct-batch-effect --batch-group-by Source -run-louvain --run-tsne example_10x.h5 example
	scrtools de_analysis --labels louvain_labels -p 20 --fisher example.h5ad example_de.xlsx
	scrtools annotate_cluster example.h5ad example.anno.txt
	scrtools plot composition --cluster-labels louvain_labels --attribute Donor --style normalized --not-stacked example.h5ad example.composition.png
	scrtools plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad example.scatter.png
	scrtools iplot --attribute louvain_labels diffmap_pca example.h5ad example.diffmap.html

The above analysis will give you tSNE, louvain cluster labels and diffusion maps in ``example.h5ad``. You can investigate donor-specific effects by looking at ``example.composition.png``. ``example.scatter.png`` plotted tSNE colored by louvain_labels and Donor info side-by-side. You can explore the diffusion map in 3D by looking at ``example.diffmap.html``. This html maps all diffusion components into 3D using PCA.

If you want to perform subcluster analysis by combining cluster 1 and 3, run the following command::

	scrtools subcluster -p 20 --correct-batch-effect example.h5ad 1,3 example_sub


---------------------------------


``scrtools aggregate_matrix``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step for single cell analysis is to generate one count matrix from cellranger's channel-specific count matrices. ``scrtools aggregate_matrix`` allows aggregating arbitrary matrices with the help of a *CSV* file.

Type::

	scrtools aggregate_matrix -h

to see the usage information::

	Usage:
  		scrtools aggregate_matrix <csv_file> <output_name> [--genome <genome> --restriction <restriction>... --attributes <attributes> --google-cloud]
  		scrtools aggregate_matrix -h

* Arguments:

	csv_file
		Input csv-formatted file containing information of each 10x channel. Each row must contain at least 3 columns --- Sample, sample name; Location, location of the channel-specific count matrix in 10x format (e.g. /sample/filtered_gene_bc_matrices_h5.h5); Reference, genome reference used for 10x cellranger. See below for an example csv::

			Sample,Source,Platform,Donor,Reference,Location
 			sample_1,bone_marrow,NextSeq,1,GRCh38,/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5
			sample_2,bone_marrow,NextSeq,2,GRCh38,/my_dir/sample_2/filtered_gene_bc_matrices_h5.h5
			sample_3,pbmc,NextSeq,1,GRCh38,/my_dir/sample_3/filtered_gene_bc_matrices_h5.h5
			sample_4,pbmc,NextSeq,2,GRCh38,/my_dir/sample_4/filtered_gene_bc_matrices_h5.h5

	output_name
		The output file name.

* Options:
	
	-\\-genome <genome>
		Genome reference. [default: GRCh38]

	-\\-restriction <restriction>...
		Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value or name:~value,..,value, where ~ refers to not. You can specifiy multiple restrictions by setting this option multiple times.

	-\\-attributes <attributes>
		Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.

	-\\-google-cloud
		If files are stored in google cloud. Assuming google cloud sdk is installed.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_name_10x.h5
		A 10x-formatted HDF5 file containing the count matrix and associated attributes.

* Examples::

	scrtools aggregate_matrix --genome GRCh38 --restriction Source:pbmc --restriction Donor:1 --attributes Source,Platform,Donor example.csv example


---------------------------------


``scrtools cluster``
^^^^^^^^^^^^^^^^^^^^

Once we collected the count matrix ``example_10x.h5``, we can perform single cell analysis using ``scrtools cluster``.

Type::

	scrtools cluster -h

to see the usage information::

	Usage:
		scrtools cluster [options] <input_file> <output_name>
		scrtools cluster -h

* Arguments:

	input_file
		Input file in 10x format. If first-pass analysis has been performed, but you want to run some additional analysis, you could also pass a h5ad-formatted file.

	output_name      
		Output file name. All outputs will use it as the prefix.

* Options:

	\-p <number>, -\\-threads <number>
		Number of threads. [default: 1]

	-\\-genome <genome>
		Genome name. [default: GRCh38]

	-\\-processed
		Input file is processed and thus no PCA & diffmap will be run.

  	-\\-output-filtration-results <spreadsheet>
		Output filtration results into <spreadsheet>.

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

	-\\-mito-prefix <prefix>
		Prefix for mitochondrial genes. [default: MT-]

	-\\-percent-mito <ratio>
		Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]

	-\\-gene-percent-cells <ratio>
		Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]

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

	-\\-run-kmeans
		Run KMeans clustering algorithm on diffusion components.

	-\\-kmeans-n-clusters <number>
		Target at <number> clusters for K means. [default: 20]

	-\\-run-hdbscan
		Run hdbscan clustering algorithm on diffusion components.

	-\\-hdbscan-min-cluster-size <number>
		Minimum cluster size for hdbscan. [default: 50]

	-\\-hdbscan-min-samples <number>
		Minimum number of samples for hdbscan. [default: 50]

	-\\-run-approximated-louvain
		Run approximated louvain clustering algorithm.

	-\\-approx-louvain-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-approx-louvain-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-approx-louvain-resolution <resolution>.
		Resolution parameter for louvain. [default: 1.3]

	-\\-run-tsne
		Run multi-core tSNE for visualization.

	-\\-tsne-perplexity <perplexity>
		tSNE's perplexity parameter. [default: 30]

  	-\\-run-fitsne
  		Run FItSNE for visualization.

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
		Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, tSNE and diffusion map coordinates are stored in the 'obsm' field.

	output_name.loom
		Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.

* Examples::

	scrtools cluster -p 20 --correct-batch-effect --run-louvain --run-tsne example_10x.h5 example



---------------------------------


``scrtools de_analysis``
^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the clusters, we can detect markers using ``scrtools de_analysis``.

Type::

	scrtools de_analysis -h

to see the usage information::

	Usage:
		scrtools de_analysis [--labels <attr> -p <threads> --alpha <alpha> --fisher --mwu --roc] <input_h5ad_file> <output_spreadsheet>
		scrtools de_analysis -h

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

	scrtools de_analysis --labels louvain_labels -p 20 --fisher --mwu --roc example.h5ad example_de.xlsx


---------------------------------


``scrtools annotate_cluster``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once we have the DE results, we could optionally identify putative cell types for each cluster using ``scrtools annotate_cluster``. Currently, this subcommand only works for human and mouse immune cells.

Type::

	scrtools annotate_cluster -h

to see the usage information::

	Usage:
		scrtools annotate_cluster [--json-file <file> --minimum-report-score <score> --do-not-use-non-de-genes] <input_h5ad_file> <output_file>
		scrtools annotate_cluster -h

* Arguments:

	input_h5ad_file
		Single cell data with DE analysis done by ``scrtools de_analysis``.

	output_file
		Output annotation file.

* Options:

	-\\-json-file <file>
		JSON file for markers. Could also be ``human``/``mouse``. [default: human]

	-\\-minimum-report-score <score>
		Minimum cell type score to report a potential cell type. [default: 0.5]

	-\\-do-not-use-non-de-genes
		Do not count non DE genes as down-regulated.

	\-h, -\\-help
		Print out help information.

* Outputs:

	output_file
		This is a text file. For each cluster, all its putative cell types are listed in descending order of the cell type score. For each putative cell type, all markers support this cell type are listed. If one putative cell type has cell subtypes, all subtypes will be listed under this cell type.

* Examples::

	scrtools annotate_cluster example.h5ad example.anno.txt


---------------------------------



``scrtools plot``
^^^^^^^^^^^^^^^^^

We can make a variety of figures using ``scrtools plot``.

Type::

	scrtools plot -h

to see the usage information::

	Usage:
  		scrtools plot [options] [--restriction <restriction>...] <plot_type> <input_h5ad_file> <output_file>
		scrtools plot -h

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
		Basis for 2D plotting, chosen from 'tsne', 'fitsne', 'umap', 'pca', 'rpca', 'fle', or 'diffmap_pca'. This option is used in 'scatter', 'scatter_groups', 'scatter_genes', and 'scatter_gene_groups'. [default: tsne]

	-\\-attributes <attrs>
		<attrs> is a comma-separated list of attributes to color the basis. This option is only used in 'scatter'.

	-\\-restriction <restriction>...
		Multiple <restriction> strings for different attributes. Each <restriction> takes the format of 'attr:value,value'. Only used for scatter.

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
		Number of rows in the figure. If not set, scrtools will figure it out automatically.

	-\\-ncols <ncols>
		Number of columns in the figure. If not set, scrtools will figure it out automatically.

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

	scrtools plot composition --cluster-labels louvain_labels --attribute Donor --style normalized --not-stacked example.h5ad example.composition.png
	scrtools plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad example.scatter.png
	scrtools plot scatter_groups --cluster-labels louvain_labels --group Donor example.h5ad example.scatter_groups.png
	scrtools plot scatter_genes --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP example.h5ad example.genes.png
	scrtools plot scatter_gene_groups --gene CD8A --group Donor example.h5ad example.gene_groups.png
	scrtools plot heatmap --cluster-labels louvain_labels --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP --heatmap-title 'markers' example.h5ad example.heatmap.png


---------------------------------



``scrtools iplot``
^^^^^^^^^^^^^^^^^^

We can also make interactive plots in html format using ``scrtools iplot``. These interactive plots are very helpful if you want to explore the diffusion maps.

Type::

	scrtools iplot -h

to see the usage information::

	Usage:
		scrtools iplot --attribute <attr> [options] <basis> <input_h5ad_file> <output_html_file>
		scrtools iplot -h

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

	scrtools iplot --attribute louvain_labels tsne example.h5ad example.tsne.html
	scrtools iplot --attribute louvain_labels diffmap_pca example.h5ad example.diffmap.html


---------------------------------

``scrtools view``
^^^^^^^^^^^^^^^^^

We may want to further perform sub-cluster analysis on a subset of cells. This sub-command helps us to define the subset.

Type::

	scrtools view -h

to see the usage information::

	Usage:
		scrtools view [--show-attributes --show-gene-attributes --show-values-for-attributes <attributes>] <input_h5ad_file>
		scrtools view -h

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

	scrtools view --show-attributes example.h5ad
	scrtools view --show-gene-attributes example.h5ad
	scrtools view --show-values-for-attributes louvain_labels,Donor example.h5ad


---------------------------------


``scrtools subcluster``
^^^^^^^^^^^^^^^^^^^^^^^

If there is a subset of cells that we want to further cluster, we can run ``scrtools subcluster``. This sub-command will outputs a new h5ad file that you can run ``de_analysis``, ``plot`` and ``iplot`` on.

Type::

	scrtools subcluster -h

to see the usage information::

	Usage:
		scrtools subcluster [options] --subset-selection <subset-selection>... <input_file> <output_name>
		scrtools subcluster -h

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

	-\\-run-kmeans
		Run KMeans clustering algorithm on diffusion components.

	-\\-kmeans-n-clusters <number>
		Target at <number> clusters for K means. [default: 20]

	-\\-run-hdbscan
		Run hdbscan clustering algorithm on diffusion components.

	-\\-hdbscan-min-cluster-size <number>
		Minimum cluster size for hdbscan. [default: 50]

	-\\-hdbscan-min-samples <number>
		Minimum number of samples for hdbscan. [default: 50]

	-\\-run-approximated-louvain
		Run approximated louvain clustering algorithm.

	-\\-approx-louvain-ninit <number>
		Number of Kmeans tries. [default: 20]

	-\\-approx-louvain-nclusters <number>
		Number of clusters for Kmeans initialization. [default: 30]

	-\\-approx-louvain-resolution <resolution>.
		Resolution parameter for louvain. [default: 1.3]

	-\\-run-tsne
		Run multi-core tSNE for visualization.

	-\\-tsne-perplexity <perplexity>
		tSNE's perplexity parameter. [default: 30]

  	-\\-run-fitsne
  		Run FItSNE for visualization.

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
		Output file in h5ad format. The clustering results are stored in the 'obs' field (e.g. 'louvain_labels' for louvain cluster labels). The PCA, tSNE and diffusion map coordinates are stored in the 'obsm' field.

	output_name.loom
		Optional output. Only exists if '--output-loom' is set. output_name.h5ad in loom format for visualization.

* Examples::

	scrtools subcluster --subset_selection louvain_labels:1,3  --subset_selection Donor:1 -p 20 --correct-batch-effect example.h5ad example_sub
