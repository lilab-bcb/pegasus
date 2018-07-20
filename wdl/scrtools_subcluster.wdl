import "https://api.firecloud.org/ga4gh/v1/tools/scrtools:tasks/versions/6/plain-WDL/descriptor" as tasks
# import "../scrtools_tasks.wdl" as tasks

workflow scrtools_subcluster {
	File input_h5ad
	# google bucket, subdirectory name and results name prefix
	String output_name
	# Specify which cells will be included in the subcluster analysis. This field contains one or more <subset_selection> strings separated by ';'. Each <subset_selection> string takes the format of ‘attr:value,…,value’, which means select cells with attr in the values. If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings.
	String subset_selections
	
	# Number of cpus per scrtools job
	Int? num_cpu = 64
	# Memory size in GB
	Int? memory = 200
	# Total disk space
	Int? disk_space = 100
	# Number of preemptible tries 
	Int? preemptible = 2


	String out_name = basename(output_name)


	# for subcluster

	# If correct batch effects [default: false]
	Boolean? correct_batch_effect
	# If output loom-formatted file [default: false]
	Boolean? output_loom
	# Random number generator seed. [default: 0]
	Int? random_state
	# Number of PCs. [default: 50]
	Int? nPC
	# Number of diffusion components. [default: 50]
	Int? nDC
	# Power parameter for diffusion-based pseudotime. [default: 0.5]
	Float? diffmap_alpha
	# Number of neighbors used for constructing affinity matrix. [default: 100]
	Float? diffmap_K
	# Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.
	String? calculate_pseudotime
	# Run louvain clustering algorithm.
	Boolean? run_louvain
	# Resolution parameter for the louvain clustering algorithm. [default: 1.3]
	Float? louvain_resolution
	# Run KMeans clustering algorithm on diffusion components.
	Boolean? run_kmeans
	# Target at <number> clusters for K means. [default: 20]
	Int? kmeans_n_clusters
	# Run hdbscan clustering algorithm on diffusion components.
	Boolean? run_hdbscan
	# Minimum cluster size for hdbscan. [default: 50]
	Int? hdbscan_min_cluster_size
	# Minimum number of samples for hdbscan. [default: 50]
	Int? hdbscan_min_samples
	# Run approximated louvain clustering algorithm.
	Boolean? run_approximated_louvain
	# Number of Kmeans tries. [default: 20]
	Int? approx_louvain_ninit
	# Number of clusters for Kmeans initialization. [default: 30]
	Int? approx_louvain_nclusters
	# Resolution parameter for louvain. [default: 1.3]
	Float? approx_louvain_resolution
	# Run multi-core tSNE for visualization.
	Boolean? run_tsne
	# tSNE’s perplexity parameter. [default: 30]
	Float? tsne_perplexity
	# Run FItSNE for visualization.
	Boolean? run_fitsne
	# Run umap for visualization.
	Boolean? run_umap
	# Run umap on diffusion components.
	Boolean? umap_on_diffmap
	# K neighbors for umap. [default: 15]
	Int? umap_K
	# Umap parameter. [default: 0.1]
	Float? umap_min_dist
	# Umap parameter. [default: 1.0]
	Float? umap_spread
	# Run force-directed layout embedding.
	Boolean? run_fle
	# K neighbors for building graph for FLE. [default: 50]
	Int? fle_K
	# Number of iterations for FLE. [default: 10000]
	Int? fle_n_steps


	# for de_analysis and annotate_cluster

	# If perform de analysis
	Boolean perform_de_analysis = true
	# Specify the cluster labels used for differential expression analysis. [default: louvain_labels]
	String? cluster_labels
	# Control false discovery rate at <alpha>. [default: 0.05]
	Float? alpha
	# Calculate Fisher’s exact test.
	Boolean? fisher
	# Calculate Mann-Whitney U test.
	Boolean? mwu
	# Calculate area under cuver in ROC curve.
	Boolean? roc

	# If also annotate cell types for clusters based on DE results.
	Boolean? annotate_cluster
	# Organism, could either be "human" or "mouse" [default: human]
	String? organism
	# Minimum cell type score to report a potential cell type. [default: 0.5]
	Float? minimum_report_score


	# for plot

	# Takes the format of "label:attr,label:attr,...,label:attr". If non-empty, generate composition plot for each "label:attr" pair. "label" refers to cluster labels and "attr" refers to sample conditions.
	String? plot_composition
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored tSNEs side by side. 
	String? plot_tsne
	# Takes the format of "attr,attr,...,attr". If non-empty, generate attr colored 3D interactive plot. The 3 coordinates are the first 3 PCs of all diffusion components.
	String? plot_diffmap


	# for scp_output

	# If generate outputs required by single cell portal
	Boolean generate_scp_outputs = false

	Boolean output_dense = false


	call tasks.run_scrtools_subcluster as subcluster {
		input:
			input_h5ad = input_h5ad,
			output_name = out_name,
			subset_selections = subset_selections,
			correct_batch_effect = correct_batch_effect,
			output_loom = output_loom,
			random_state = random_state,
			nPC = nPC,
			nDC = nDC,
			diffmap_alpha = diffmap_alpha,
			diffmap_K = diffmap_K,
			calculate_pseudotime = calculate_pseudotime,
			run_louvain = run_louvain,
			louvain_resolution = louvain_resolution,
			run_kmeans = run_kmeans,
			kmeans_n_clusters = kmeans_n_clusters,
			run_hdbscan = run_hdbscan,
			hdbscan_min_cluster_size = hdbscan_min_cluster_size,
			hdbscan_min_samples = hdbscan_min_samples,
			run_approximated_louvain = run_approximated_louvain,
			approx_louvain_ninit = approx_louvain_ninit,
			approx_louvain_nclusters = approx_louvain_nclusters,
			approx_louvain_resolution = approx_louvain_resolution,
			run_tsne = run_tsne,
			tsne_perplexity = tsne_perplexity,
			run_fitsne = run_fitsne,
			run_umap = run_umap,
			umap_on_diffmap = umap_on_diffmap,
			umap_K = umap_K,
			umap_min_dist = umap_min_dist,
			umap_spread = umap_spread,
			run_fle = run_fle,
			fle_K = fle_K,
			fle_n_steps = fle_n_steps,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible
	}

	if (perform_de_analysis) {
		call tasks.run_scrtools_de_analysis as de_analysis {
			input:
				input_h5ad = subcluster.output_h5ad,
				output_name = out_name,
				labels = cluster_labels,
				alpha = alpha,
				fisher = fisher,
				mwu = mwu,
				roc = roc,
				annotate_cluster = annotate_cluster,
				organism = organism,
				minimum_report_score = minimum_report_score,
				num_cpu = num_cpu,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible
		}
	}

	if (defined(plot_composition) || defined(plot_tsne) || defined(plot_diffmap)) {
		call tasks.run_scrtools_plot as plot {
			input:
				input_h5ad = subcluster.output_h5ad,
				output_name = out_name,
				plot_composition = plot_composition,
				plot_tsne = plot_tsne,
				plot_diffmap = plot_diffmap,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible
		}	
	}

	if (generate_scp_outputs) {
		call tasks.run_scrtools_scp_output as scp_output {
			input:
				input_h5ad = subcluster.output_h5ad,
				output_name = out_name,
				output_dense = output_dense,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible				
		}
	}

	call tasks.organize_results {
		input:
			output_name = output_name,
			output_h5ad = subcluster.output_h5ad,
			output_loom_file = subcluster.output_loom_file,
			output_de_h5ad = de_analysis.output_de_h5ad,
			output_de_xlsx = de_analysis.output_de_xlsx,
			output_anno_file = de_analysis.output_anno_file,
			output_pngs = plot.output_pngs,
			output_htmls = plot.output_htmls,
			output_scp_files = scp_output.output_scp_files,
			disk_space = disk_space,
			preemptible = preemptible
	}
}
