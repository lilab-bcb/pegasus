workflow scrtools_subcluster {
	# google bucket and subdirectory name
	String data_folder
	String input_name
	String output_name
	String cluster_ids

	Int num_cpu = 64
	Boolean? correct_batch_effect
	Boolean? output_loom
	# Color tSNE by <attribute> and put it at the right side of louvain clusters.
	String? plot_by_side
	# Put legends on the tSNE.
	Boolean? legend_on_data
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

	call run_scrtools_subcluster {
		input:
			data_folder = data_folder,
			input_name = input_name,
			output_name = output_name,
			cluster_ids = cluster_ids,
			num_cpu = num_cpu,
			correct_batch_effect = correct_batch_effect,
			output_loom = output_loom,
			plot_by_side = plot_by_side,
			legend_on_data = legend_on_data,
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

task run_scrtools_subcluster {
	String data_folder
	String input_name
	String output_name
	String cluster_ids
	Int num_cpu
	Boolean? correct_batch_effect
	Boolean? output_loom
	String? plot_by_side
	Boolean? legend_on_data
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
		cluster_ids = '${cluster_ids}'
		plot_by_side = '${plot_by_side}'
		louvain_resolution = '${louvain_resolution}'		

		call_args = ['gsutil', '-m', 'cp', '${data_folder}/' + input_name + '.h5ad', '.']
		print(' '.join(call_args))
		check_call(call_args)

		call_args = ['scrtools', 'subcluster', input_name + '.h5ad', cluster_ids, output_name, '-p', '${num_cpu}']
		if '${correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if plot_by_side is not '':
			call_args.extend(['--plot-by-side', plot_by_side])
		if '${legend_on_data}' is 'true':
			call_args.append('--legend-on-data')
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
		docker: "regevlab/scrtools"
		memory: "200 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: "${num_cpu}"
		preemptible: 2
	}
}
