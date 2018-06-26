Use ``scrtools`` as a command line tool
---------------------------------------

``scrtools`` can be used as a command line tool. Type::

	scrtools -h

to see the help information::

	Usage:
		scrtools <command> [<args>...]
		scrtools -h | --help
		scrtools -v | --version

``scrtools`` has 7 sub-commands in 3 groups.

* Preprocessing:

	aggregate_matrix
		Aggregate cellranger-outputted channel-specific count matrices into a single count matrix. It also enables users to import metadata into the count matrix.

* Analyzing:
	
	cluster
		Perform first-pass analysis using the count matrix generated from 'aggregate_matrix'. This subcommand could perform low quality cell filtration, batch correction, variable gene selection, dimension reduction, diffusion map calculation, graph-based clustering, tSNE visualization. The final results will be written into h5ad-formatted file, which Seurat could load.
    
    subcluster
    	Perform sub-cluster analyses based on one or several clusters obtained by 'cluster'.
  		
    de_analysis
    	Detect markers for each cluster by performing differential expression analysis per cluster (within cluster vs. outside cluster). DE tests include Welch's t-test, Fisher's exact test, Mann-Whitney U test. It can also calculate AUROC values for each gene.
    
    annotate_cluster
    	This subcommand is used to automatically annotate cell types for each cluster based on existing markers. Currently, it only works for human and mouse immune cells.

* Plotting:

	plot
		Make static plots, which includes plotting tSNEs by cluster labels and different groups.
			
	iplot
		Make interactive plots using plotly. The outputs are HTML pages. You can visualize diffusion maps with this sub-command.



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
		Input csv-formatted file containing information of each 10x channel. Each row must contain at least 3 columns --- Sample, sample name; Location, folder that contains the count matrices (e.g. filtered_gene_bc_matrices_h5.h5); Reference, genome reference used for 10x cellranger. See below for an example csv::

			Sample,Source,Platform,Donor,Reference,Location
 			S1,bone_marrow,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S1
			S2,bone_marrow,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S2
			S3,pbmc,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S3
			S4,pbmc,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S4

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

	scrtools aggregate_matrix --genome GRCh38 --restriction Source:pbmc --restriction Donor:1 --attributes Source,Platform example.csv example


---------------------------------


``scrtools cluster``
^^^^^^^^^^^^^^^^^^^^

Once we collected the count matrix ``example_10x.h5``, we can perform single cell analysis using ``scrtools cluster``.

Type::

	scrtools aggregate_matrix -h

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
		Group batches according to <expression>. If <expression> is None, assume all channels are of one group.
  
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

	-\\-calculate-pseudotime <roots>
		Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a list of cell barcodes.

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


``scrtools 