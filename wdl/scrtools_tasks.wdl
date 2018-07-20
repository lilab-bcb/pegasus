workflow scrtools_tasks {
}

task run_scrtools_aggregate_matrices {
	File input_count_matrix_csv
	String output_name
	Int memory
	Int disk_space
	Int preemptible
	String? genome
	String? restrictions
	String? attributes

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		call_args = ['scrtools', 'aggregate_matrix', '${input_count_matrix_csv}', '${output_name}', '--google-cloud']
		if '${genome}' is not '':
			call_args.extend(['--genome', '${genome}'])
		if '${restrictions}' is not '':
			ress = '${restrictions}'.split(';')
			for res in ress:
				call_args.extend(['--restriction', res])
		if '${attributes}' is not '':
			call_args.extend(['--attributes', '${attributes}'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_10x_h5 = '${output_name}_10x.h5'
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_scrtools_cluster {
	File input_10x_file
	String output_name
	Int num_cpu
	Int memory
	Int disk_space
	Int preemptible
	String? genome
	Boolean? output_filtration_results
	Boolean? output_loom
	Boolean? correct_batch_effect
	String? batch_group_by
	Int? min_genes
	Int? max_genes
	String? mito_prefix
	Float? percent_mito
	Float? gene_percent_cells
	Float? counts_per_cell_after
	Int? random_state
	Boolean? run_uncentered_pca
	Boolean? no_variable_gene_selection
	Boolean? no_submat_to_dense
	Int? nPC
	Int? nDC
	Float? diffmap_alpha
	Float? diffmap_K
	Boolean? run_louvain
	Float? louvain_resolution
	String? louvain_affinity
	Boolean? run_kmeans
	Int? kmeans_n_clusters
	Boolean? run_hdbscan
	Int? hdbscan_min_cluster_size
	Int? hdbscan_min_samples
	Boolean? run_approximated_louvain
	Int? approx_louvain_ninit
	Int? approx_louvain_nclusters
	Float? approx_louvain_resolution
	Boolean? run_tsne
	Float? tsne_perplexity
	Boolean? run_fitsne
	Boolean? run_umap
	Boolean? umap_on_diffmap
	Int? umap_K
	Float? umap_min_dist
	Float? umap_spread
	Boolean? run_fle
	Int? fle_K
	Int? fle_n_steps
	String? fle_affinity

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['scrtools', 'cluster', '${input_10x_file}', '${output_name}', '-p', '${num_cpu}']
		if '${genome}' is not '':
			call_args.extend(['--genome', '${genome}'])
		if '${output_filtration_results}' is 'true':
			call_args.extend(['--output-filtration-results', '${output_name}' + '.filt.xlsx'])
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if '${correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
			if '${batch_group_by}' is not '':
				call_args.extend(['--batch-group-by', '${batch_group_by}'])
		if '${min_genes}' is not '':
			call_args.extend(['--min-genes', '${min_genes}'])
		if '${max_genes}' is not '':
			call_args.extend(['--max-genes', '${max_genes}'])
		if '${mito_prefix}' is not '':
			call_args.extend(['--mito-prefix', '${mito_prefix}'])
		if '${percent_mito}' is not '' :
			call_args.extend(['--percent-mito', '${percent_mito}'])
		if '${gene_percent_cells}' is not '':
			call_args.extend(['--gene-percent-cells', '${gene_percent_cells}'])
		if '${counts_per_cell_after}' is not '':
			call_args.extend(['--counts-per-cell-after', '${counts_per_cell_after}'])
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${run_uncentered_pca}' is 'true':
			call_args.append('--run-uncentered-pca')
		if '${no_variable_gene_selection}' is 'true':
			call_args.append('--no-variable-gene-selection')
		if '${no_submat_to_dense}' is 'true':
			call_args.append('--no-submat-to-dense')
		if '${nPC}' is not '':
			call_args.extend(['--nPC', '${nPC}'])
		if '${nDC}' is not '':
			call_args.extend(['--nDC', '${nDC}'])
		if '${diffmap_alpha}' is not '':
			call_args.extend(['--diffmap-alpha', '${diffmap_alpha}'])
		if '${diffmap_K}' is not '':
			call_args.extend(['--diffmap-K', '${diffmap_K}'])
		if '${run_louvain}' is 'true':
			call_args.append('--run-louvain')
		if '${louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '${louvain_resolution}'])
		if '${louvain_affinity}' is not '':
			call_args.extend(['--louvain-affinity', '${louvain_affinity}'])
		if '${run_kmeans}' is 'true':
			call_args.append('--run-kmeans')
		if '${kmeans_n_clusters}' is not '':
			call_args.extend(['--kmeans-n-clusters', '${kmeans_n_clusters}'])
		if '${run_hdbscan}' is 'true':
			call_args.append('--run-hdbscan')
		if '${hdbscan_min_cluster_size}' is not '':
			call_args.extend(['--hdbscan-min-cluster-size', '${hdbscan_min_cluster_size}'])
		if '${hdbscan_min_samples}' is not '':
			call_args.extend(['--hdbscan-min-samples', '${hdbscan_min_samples}'])
		if '${run_approximated_louvain}' is 'true':
			call_args.append('--run-approximated-louvain')
		if '${approx_louvain_ninit}' is not '':
			call_args.extend(['--approx-louvain-ninit', '${approx_louvain_ninit}'])
		if '${approx_louvain_nclusters}' is not '':
			call_args.extend(['--approx-louvain-nclusters', '${approx_louvain_nclusters}'])
		if '${approx_louvain_resolution}' is not '':
			call_args.extend(['--approx-louvain-resolution', '${approx_louvain_resolution}'])
		if '${run_tsne}' is 'true':
			call_args.append('--run-tsne')
		if '${tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '${tsne_perplexity}'])
		if '${run_fitsne}' is 'true':
			call_args.append('--run-fitsne')
		if '${run_umap}' is 'true':
			call_args.append('--run-umap')
		if '${umap_on_diffmap}' is 'true':
			call_args.append('--umap-on-diffmap')
		if '${umap_K}' is not '':
			call_args.extend(['--umap-K', '${umap_K}'])
		if '${umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '${umap_min_dist}'])
		if '${umap_spread}' is not '':
			call_args.extend(['--umap-spread', '${umap_spread}'])
		if '${run_fle}' is 'true':
			call_args.append('--run-fle')
		if '${fle_K}' is not '':
			call_args.extend(['--fle-K', '${fle_K}'])
		if '${fle_n_steps}' is not '':
			call_args.extend(['--fle-n-steps', '${fle_n_steps}'])
		if '${fle_affinity}' is not '':
			call_args.extend(['--fle-affinity', '${fle_affinity}'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_h5ad = "${output_name}.h5ad"
		Array[File] output_filt_xlsx = glob("${output_name}.filt.xlsx")
		Array[File] output_loom_file = glob("${output_name}.loom")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_scrtools_de_analysis {
	File input_h5ad
	String output_name
	Int num_cpu
	Int memory
	Int disk_space
	Int preemptible
	String? labels
	Float? alpha
	Boolean? fisher
	Boolean? mwu
	Boolean? roc

	Boolean? annotate_cluster
	String? organism
	Float? minimum_report_score

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['mv', '-f', '${input_h5ad}', '${output_name}.h5ad']
		print(' '.join(call_args))
		check_call(call_args)			
		call_args = ['scrtools', 'de_analysis', '${output_name}.h5ad', '${output_name}.de.xlsx', '-p', '${num_cpu}']
		if '${labels}' is not '':
			call_args.extend(['--labels', '${labels}'])
		if '${alpha}' is not '':
			call_args.extend(['--alpha', '${alpha}'])
		if '${fisher}' is 'true':
			call_args.append('--fisher')
		if '${mwu}' is 'true':
			call_args.append('--mwu')
		if '${roc}' is 'true':
			call_args.append('--roc')
		print(' '.join(call_args))
		check_call(call_args)
		if '${annotate_cluster}' is 'true':
			call_args = ['scrtools', 'annotate_cluster', '${output_name}.h5ad', '${output_name}' + '.anno.txt']
			if '${organism}' is not '':
				call_args.extend(['--json-file', '${organism}'])
			if '${minimum_report_score}' is not '':
				call_args.extend(['--minimum-report-score', '${minimum_report_score}'])
			print(' '.join(call_args))
			check_call(call_args)			
		CODE
	}

	output {
		File output_de_h5ad = "${output_name}.h5ad"
		File output_de_xlsx = "${output_name}.de.xlsx"
		Array[File] output_anno_file = glob("${output_name}.anno.txt")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_scrtools_plot {
	File input_h5ad
	String output_name
	Int memory
	Int disk_space
	Int preemptible
	String? plot_composition
	String? plot_tsne
	String? plot_diffmap

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		if '${plot_composition}' is not '':
			pairs = '${plot_composition}'.split(',')
			for pair in pairs:
				lab, attr = pair.split(':')
				call_args = ['scrtools', 'plot', 'composition', '--cluster-labels', lab, '--attribute', attr, '--style', 'normalized', '--not-stacked', '${input_h5ad}', '${output_name}.' + lab + '+' + attr + '.composition.png']
				print(' '.join(call_args))
				check_call(call_args)
		if '${plot_tsne}' is not '':
			call_args = ['scrtools', 'plot', 'scatter', '--attributes', '${plot_tsne}', '${input_h5ad}', '${output_name}.tsne.png']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_diffmap}' is not '':
			attrs = '${plot_diffmap}'.split(',')
			for attr in attrs:
				call_args = ['scrtools', 'iplot', '--attribute', attr, 'diffmap_pca', '${input_h5ad}', '${output_name}.' + attr + '.diffmap_pca.html']
				print(' '.join(call_args))
				check_call(call_args)
		CODE
	}

	output {
		Array[File] output_pngs = glob("*.png")
		Array[File] output_htmls = glob("*.html")
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_scrtools_scp_output {
	File input_h5ad
	String output_name
	Boolean output_dense
	Int memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp
		scrtools scp_output ${true='--dense' false='' output_dense} ${input_h5ad} ${output_name}
	}

	output {
		Array[File] output_scp_files = glob("${output_name}.scp.*")
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_scrtools_subcluster {
	File input_h5ad
	String output_name
	Int num_cpu
	Int memory
	Int disk_space
	Int preemptible
	String subset_selections 
	Boolean? correct_batch_effect
	Boolean? output_loom
	Int? random_state
	Boolean? run_uncentered_pca
	Boolean? no_variable_gene_selection
	Boolean? no_submat_to_dense
	Int? nPC
	Int? nDC
	Float? diffmap_alpha
	Float? diffmap_K
	String? calculate_pseudotime
	Boolean? run_louvain
	Float? louvain_resolution
	String? louvain_affinity
	Boolean? run_kmeans
	Int? kmeans_n_clusters
	Boolean? run_hdbscan
	Int? hdbscan_min_cluster_size
	Int? hdbscan_min_samples
	Boolean? run_approximated_louvain
	Int? approx_louvain_ninit
	Int? approx_louvain_nclusters
	Float? approx_louvain_resolution
	Boolean? run_tsne
	Float? tsne_perplexity
	Boolean? run_fitsne
	Boolean? run_umap
	Boolean? umap_on_diffmap
	Int? umap_K
	Float? umap_min_dist
	Float? umap_spread
	Boolean? run_fle
	Int? fle_K
	Int? fle_n_steps
	String? fle_affinity

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['scrtools', 'subcluster', '${input_h5ad}', '${output_name}', '-p', '${num_cpu}']
		if '${subset_selections}' is not '':
			sels = '${subset_selections}'.split(';')
			for sel in sels:
				call_args.extend(['--subset-selection', sel])
		if '${correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${run_uncentered_pca}' is 'true':
			call_args.append('--run-uncentered-pca')
		if '${no_variable_gene_selection}' is 'true':
			call_args.append('--no-variable-gene-selection')
		if '${no_submat_to_dense}' is 'true':
			call_args.append('--no-submat-to-dense')
		if '${nPC}' is not '':
			call_args.extend(['--nPC', '${nPC}'])
		if '${nDC}' is not '':
			call_args.extend(['--nDC', '${nDC}'])
		if '${diffmap_alpha}' is not '':
			call_args.extend(['--diffmap-alpha', '${diffmap_alpha}'])
		if '${diffmap_K}' is not '':
			call_args.extend(['--diffmap-K', '${diffmap_K}'])
		if '${calculate_pseudotime}' is not '':
			call_args.extend(['--calculate-pseudotime', '${calculate_pseudotime}'])
		if '${run_louvain}' is 'true':
			call_args.append('--run-louvain')
		if '${louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '${louvain_resolution}'])
		if '${louvain_affinity}' is not '':
			call_args.extend(['--louvain-affinity', '${louvain_affinity}'])
		if '${run_kmeans}' is 'true':
			call_args.append('--run-kmeans')
		if '${kmeans_n_clusters}' is not '':
			call_args.extend(['--kmeans-n-clusters', '${kmeans_n_clusters}'])
		if '${run_hdbscan}' is 'true':
			call_args.append('--run-hdbscan')
		if '${hdbscan_min_cluster_size}' is not '':
			call_args.extend(['--hdbscan-min-cluster-size', '${hdbscan_min_cluster_size}'])
		if '${hdbscan_min_samples}' is not '':
			call_args.extend(['--hdbscan-min-samples', '${hdbscan_min_samples}'])
		if '${run_approximated_louvain}' is 'true':
			call_args.append('--run-approximated-louvain')
		if '${approx_louvain_ninit}' is not '':
			call_args.extend(['--approx-louvain-ninit', '${approx_louvain_ninit}'])
		if '${approx_louvain_nclusters}' is not '':
			call_args.extend(['--approx-louvain-nclusters', '${approx_louvain_nclusters}'])
		if '${approx_louvain_resolution}' is not '':
			call_args.extend(['--approx-louvain-resolution', '${approx_louvain_resolution}'])
		if '${run_tsne}' is 'true':
			call_args.append('--run-tsne')
		if '${tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '${tsne_perplexity}'])
		if '${run_fitsne}' is 'true':
			call_args.append('--run-fitsne')
		if '${run_umap}' is 'true':
			call_args.append('--run-umap')
		if '${umap_on_diffmap}' is 'true':
			call_args.append('--umap-on-diffmap')
		if '${umap_K}' is not '':
			call_args.extend(['--umap-K', '${umap_K}'])
		if '${umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '${umap_min_dist}'])
		if '${umap_spread}' is not '':
			call_args.extend(['--umap-spread', '${umap_spread}'])
		if '${run_fle}' is 'true':
			call_args.append('--run-fle')
		if '${fle_K}' is not '':
			call_args.extend(['--fle-K', '${fle_K}'])
		if '${fle_n_steps}' is not '':
			call_args.extend(['--fle-n-steps', '${fle_n_steps}'])
		if '${fle_affinity}' is not '':
			call_args.extend(['--fle-affinity', '${fle_affinity}'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_h5ad = "${output_name}.h5ad"
		Array[File] output_loom_file = glob("${output_name}.loom")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task organize_results {
	String output_name
	Int disk_space
	Int preemptible
	File? output_10x_h5
	File? output_h5ad
	Array[File]? output_filt_xlsx
	Array[File]? output_loom_file
	File? output_de_h5ad
	File? output_de_xlsx
	Array[File]? output_anno_file
	Array[File]? output_pngs
	Array[File]? output_htmls
	Array[File]? output_scp_files

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import os
		from subprocess import check_call
		dest = os.path.dirname('${output_name}') + '/'
		files = ['${output_10x_h5}', '${sep=" " output_filt_xlsx}', '${sep=" " output_loom_file}', '${output_de_xlsx}', '${sep=" " output_anno_file}']
		files.append('${output_h5ad}' if '${output_de_h5ad}' is '' else '${output_de_h5ad}')
		files.extend('${sep="," output_pngs}'.split(','))
		files.extend('${sep="," output_htmls}'.split(','))
		files.extend('${sep="," output_scp_files}'.split(','))
		for file in files:
			if file is not '':
				# call_args = ['cp', file, dest]
				call_args = ['gsutil', '-q', 'cp', file, dest]
				print(' '.join(call_args))
				check_call(call_args)
		CODE
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "30 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}
