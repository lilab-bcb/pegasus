workflow scrtools_aggregate_matrices {
	File input_count_matrix_csv
	# google bucket and subdirectory name
	String output_folder
	# output file name prefix
	String output_name
	
	String? genome
	# Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. Multiple restrictions are separated by ';'
	String? restrictions
	# Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.
	String? attributes

	# When perform batch correction, perform it in each group separately, should be in the format of 'value' or 'attribute=value' or 'attr1+attr2+....'
	String? groupby

	Int? diskSpace = 100

	call run_scrtools_aggregate_matrices {
		input:
			input_count_matrix_csv = input_count_matrix_csv,
			output_folder = output_folder,
			output_name = output_name,
			genome = genome,
			restrictions = restrictions,
			attributes = attributes,
			groupby = groupby,
			diskSpace = diskSpace
	}
}

task run_scrtools_aggregate_matrices {
	File input_count_matrix_csv
	String output_folder
	String output_name
	String? genome
	String? restrictions
	String? attributes
	String? groupby
	Int? diskSpace

	command {
		set -e
		export TMPDIR=/tmp
		python <<CODE
		from subprocess import check_call

		genome = '${genome}'
		restrictions = '${restrictions}'
		attributes = '${attributes}'

		call_args = ['scrtools', 'aggregate_matrix', '${input_count_matrix_csv}', '${output_name}', '--google-cloud']
		if genome is not '':
			call_args.extend(['--genome', genome])
		if restrictions is not '':
			ress = restrictions.split(';')
			for res in ress:
				call_args.extend(['--restriction', res])
		if attributes is not '':
			call_args.extend(['--attributes', attributes])

		print(' '.join(call_args))
		check_call(call_args)

		if '${groupby}' is not '':
			call_args = ['scrtools', 'add_attribute', '${output_name}.attr.csv', 'GroupBy:${groupby}']
			print(' '.join(call_args))
			check_call(call_args)
		CODE
		gsutil -m mv ${output_name}_10x.h5 ${output_folder}
		gsutil -m mv ${output_name}.attr.csv ${output_folder}
	}

	output {
		String output_10x_hdf = '${output_folder}/${output_name}_10x.h5'
		String output_attr_file = '${output_folder}/${output_name}.attr.csv'
	}

	runtime {
		docker: "regevlab/scrtools"
		memory: "30 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: 1
		preemptible: 2
	}
}
