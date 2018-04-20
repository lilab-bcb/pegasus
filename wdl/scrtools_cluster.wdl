workflow scrtools_cluster {
	# google bucket and subdirectory name
	String data_folder
	String input_name
	String output_name
	Int num_cpu = 64
	
	String? genome
	
	Boolean? output_filtration_results
	Boolean? output_loom
	Boolean? correct_batch_effect
	# Group batches by. If this option is off, make all batches as a single group.
	String? batch_group_by

	# Import attributes contained in the comma-separated list into the analysis object.
	String? import_attributes
	# Color tSNE by <attribute> and put it at the right side of louvain clusters.
	String? plot_by_side
	# Put legends on the tSNE.
	Boolean? legend_on_data

	# Only keep cells with at least <number> of genes. [default: 500]
	Int? min_genes
	# Only keep cells with less than <number> of genes. [default: 6000]
	Int? max_genes
	# Prefix for mitochondrial genes. [default: MT-]
	String? mito_prefix
	# Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]
	Float? percent_mito

	# Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]
	Float? gene_percent_cells

	# Total counts per cell after normalization. [default: 1e5]
	Float? counts_per_cell_after
	# Resolution parameter for the louvain clustering algorithm. [default: 1.3]
	Float? louvain_resolution

	Boolean? de_analysis = true
	# Minimum fold change in either percentage (fisher test) or log expression (t test) to report a DE gene. [default: 1.5]
	Float? fold_change
	# <attr> used as cluster labels. [default: louvain_labels]
	String? labels

	Boolean? plot_composition
	String? figure_size
	Boolean? plot_diffusion_map

	Int? diskSpace = 250

	call run_scrtools_cluster {
		input:
			data_folder = data_folder,
			input_name = input_name,
			output_name = output_name,
			num_cpu = num_cpu,
			genome = genome,
			output_filtration_results = output_filtration_results,
			output_loom = output_loom,
			correct_batch_effect = correct_batch_effect,
			batch_group_by = batch_group_by,
			import_attributes = import_attributes,
			plot_by_side = plot_by_side,
			legend_on_data = legend_on_data,
			min_genes = min_genes,
			max_genes = max_genes,
			mito_prefix = mito_prefix,
			percent_mito = percent_mito,
			gene_percent_cells = gene_percent_cells,
			counts_per_cell_after = counts_per_cell_after,
			louvain_resolution = louvain_resolution,
			de_analysis = de_analysis,
			fold_change = fold_change,
			labels = labels,
			plot_composition = plot_composition,
			figure_size = figure_size,
			plot_diffusion_map = plot_diffusion_map,
			diskSpace = diskSpace
	}
}

task run_scrtools_cluster {
	String data_folder
	String input_name
	String output_name
	Int num_cpu
	String? genome
	Boolean? output_filtration_results
	Boolean? output_loom
	Boolean? correct_batch_effect
	String? batch_group_by
	String? import_attributes
	String? plot_by_side
	Boolean? legend_on_data
	Int? min_genes
	Int? max_genes
	String? mito_prefix
	Float? percent_mito
	Float? gene_percent_cells
	Float? counts_per_cell_after
	Float? louvain_resolution
	Boolean? de_analysis
	Float? fold_change
	String? labels
	Boolean? plot_composition
	String? figure_size
	Boolean? plot_diffusion_map	
	Int? diskSpace

	command {
		set -e
		export TMPDIR=/tmp
		python <<CODE
		from subprocess import check_call

		input_name = '${input_name}'
		output_name = '${output_name}'
		genome = '${genome}'
		correct_batch = '${correct_batch_effect}' is 'true'
		batch_group_by = '${batch_group_by}'
		import_attributes = '${import_attributes}'
		plot_by_side = '${plot_by_side}'
		min_genes = '${min_genes}'
		max_genes = '${max_genes}'
		mito_prefix = '${mito_prefix}'
		percent_mito = '${percent_mito}'
		gene_percent_cells = '${gene_percent_cells}'
		counts_per_cell_after = '${counts_per_cell_after}'
		louvain_resolution = '${louvain_resolution}'		

		attrs = set()
		if import_attributes is not '':
			attrs = attrs | set(import_attributes.split(','))
		if batch_group_by is not '':
			attrs = attrs | set(batch_group_by.split(','))
		if plot_by_side is not '':
			attrs = attrs | set(plot_by_side.split(','))
		import_attributes = ','.join(attrs)

		call_args = ['gsutil', '-m', 'cp', '${data_folder}/' + input_name + '_10x.h5', '.']
		print(' '.join(call_args))
		check_call(call_args)
		call_args = ['gsutil', '-m', 'cp', '${data_folder}/' + input_name + '.attr.csv', '.']
		print(' '.join(call_args))
		check_call(call_args)

		call_args = ['scrtools', 'cluster', input_name, output_name, '-p', '${num_cpu}']
		if genome is not '':
			call_args.extend(['--genome', genome])
		if '${output_filtration_results}' is 'true':
			call_args.extend(['--output-filtration-results', output_name + '.filt.xlsx'])
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if correct_batch:
			call_args.append('--correct-batch-effect')
		if correct_batch and (batch_group_by is not ''):
			call_args.extend(['--batch-group-by', batch_group_by])
		if import_attributes is not '':
			call_args.extend(['--import-attributes', import_attributes])
		if plot_by_side is not '':
			call_args.extend(['--plot-by-side', plot_by_side])
		if '${legend_on_data}' is 'true':
			call_args.append('--legend-on-data')
		if min_genes is not '' :
			call_args.extend(['--min-genes', min_genes])
		if max_genes is not '' :
			call_args.extend(['--max-genes', max_genes])
		if mito_prefix is not '' :
			call_args.extend(['--mito-prefix', mito_prefix])
		if percent_mito is not '' :
			call_args.extend(['--percent-mito', percent_mito])
		if gene_percent_cells is not '' :
			call_args.extend(['--gene-percent-cells', gene_percent_cells])
		if counts_per_cell_after is not '' :
			call_args.extend(['--counts-per-cell-after', counts_per_cell_after])
		if louvain_resolution is not '' :
			call_args.extend(['--louvain-resolution', louvain_resolution])

		print(' '.join(call_args))
		check_call(call_args)

		if '${de_analysis}' is 'true':
			fold_change = '${fold_change}'
			labels = '${labels}'

			call_args = ['scrtools', 'de_analysis', output_name + '.h5ad', output_name]
			if fold_change is not '':
				call_args.extend(['--fold-change', fold_change])
			if labels is not '':
				call_args.extend(['--labels', labels])
			print(' '.join(call_args))
			check_call(call_args)

		if '${plot_composition}' is 'true':
			figure_size = '${figure_size}'
			call_args = ['scrtools', 'plot', 'composition', '--attribute', plot_by_side, output_name + '.h5ad', output_name + '.composition.frequency']
			if figure_size is not '':
				call_args.extend(['--sizes', figure_size])
			print(' '.join(call_args))
			check_call(call_args)
			call_args = ['scrtools', 'plot', 'composition', '--style', 'normalized', '--not-stacked', '--attribute', plot_by_side, output_name + '.h5ad', output_name + '.composition.normalized']
			if figure_size is not '':
				call_args.extend(['--sizes', figure_size])
			print(' '.join(call_args))
			check_call(call_args)

		if '${plot_diffusion_map}' is 'true':
			call_args = ['scrtools', 'iplot', 'diffmap', '--attribute', plot_by_side, output_name + '_var.h5ad', output_name + '.diffmap.html']
			print(' '.join(call_args))
			check_call(call_args)
		CODE
		gsutil -m mv ${output_name}.h5ad ${data_folder}
		gsutil -m mv ${output_name}_var.h5ad ${data_folder}
		gsutil -m mv tsne.${output_name}.png ${data_folder}/${output_name}.tsne.png
		if [ -f ${output_name}.filt.xlsx ]
		then
			gsutil -m mv ${output_name}.filt.xlsx ${data_folder}
		fi
		if [ -f ${output_name}.loom ]
		then
			gsutil -m mv ${output_name}_var.loom ${data_folder}
			gsutil -m mv ${output_name}.loom ${data_folder}
		fi
		if [ -f ${output_name}_de.h5ad ]
		then
			gsutil -m mv ${output_name}_de.h5ad ${data_folder}
			gsutil -m mv ${output_name}_de_analysis_fisher.xlsx ${data_folder}
			gsutil -m mv ${output_name}_de_analysis_t.xlsx ${data_folder}
		fi
		if [ -f ${output_name}.composition.frequency.png ]
		then
			gsutil -m mv ${output_name}.composition.frequency.png ${data_folder}
			gsutil -m mv ${output_name}.composition.normalized.png ${data_folder}
		fi
		if [ -f ${output_name}.diffmap.html ]
		then
			gsutil -m mv ${output_name}.diffmap.html ${data_folder}
		fi
	}

	output {
		String output_h5ad = '${data_folder}/${output_name}.h5ad'
		String output_var_h5ad = '${data_folder}/${output_name}_var.h5ad'
		String output_tsne_png = '${data_folder}/${output_name}.tsne.png'
		String output_filt_xlsx = '${data_folder}/${output_name}.filt.xlsx'
		String output_loom_file = '${data_folder}/${output_name}.loom'
		String output_var_loom = '${data_folder}/${output_name}_var.loom'
		String output_de_h5ad = '${data_folder}/${output_name}_de.h5ad'
		String output_de_fisher = '${data_folder}/${output_name}_de_analysis_fisher.xlsx'
		String output_de_t = '${data_folder}/${output_name}_de_analysis_t.xlsx'
		String output_composition_frequency_png = '${data_folder}/${output_name}.composition.frequency.png'
		String output_composition_normalized_png = '${data_folder}/${output_name}.composition.normalized.png'
		String output_diffmap_html = '${data_folder}/${output_name}.diffmap.html'
	}

	runtime {
		docker: "bigbadbo/scrtools"
		memory: "200 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: "${num_cpu}"
		preemptible: 2
	}
}
