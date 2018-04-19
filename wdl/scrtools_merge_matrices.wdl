workflow scrtools_merge_matrices {
	# google bucket and subdirectory name
	String data_folder
	# comma-separated 
	String input_names
	String output_name

	String? genome
	# A comma-separated list of symbols representing each input matrix
	String? symbols
	# A comma-separated list of attributes. When merging matrices, the matrix symbol will be added in front of the attributes in the list.
	String? attributes

	Int? diskSpace = 100

	call run_scrtools_merge_matrices {
		input:
			data_folder = data_folder,
			input_names = input_names,
			output_name = output_name,
			genome = genome,
			symbols = symbols,
			attributes = attributes,
			diskSpace = diskSpace
	}
}

task run_scrtools_merge_matrices {
	String data_folder
	String input_names
	String output_name
	String? genome
	String? symbols
	String? attributes
	Int? diskSpace

	command {
		set -e
		export TMPDIR=/tmp
		python <<CODE
		from subprocess import check_call

		input_names = '${input_names}'.split(',')
		for input_name in input_names:
			call_args = ['gsutil', '-m', 'cp', '${data_folder}/' + input_name + '_10x.h5', '.']
			check_call(call_args)
			call_args = ['gsutil', '-m', 'cp', '${data_folder}/' + input_name + '.attr.csv', '.']
			check_call(call_args)

		genome = '${genome}'
		symbols = '${symbols}'
		attributes = '${attributes}'

		call_args = ['scrtools', 'merge_matrix']
		for input_name in input_names:
			call_args.extend(['-i', input_name])
		call_args.append('${output_name}')

		if genome is not '':
			call_args.extend(['--genome', genome])
		if symbols is not '':
			call_args.extend(['--symbols', symbols])
		if attributes is not '':
			call_args.extend(['--attributes', attributes])

		print(' '.join(call_args))
		check_call(call_args)
		CODE
		gsutil -m mv ${output_name}* ${data_folder}
	}

	output {
		String output_10x_hdf = '${data_folder}/${output_name}_10x.h5'
		String output_attr_file = '${data_folder}/${output_name}.attr.csv'
	}

	runtime {
		docker: "bigbadbo/scrtools"
		memory: "30 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: 1
		preemptible: 2
	}
}
