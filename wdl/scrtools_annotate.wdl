workflow scrtools_annotate {
	# google bucket and subdirectory name
	String data_folder
	# file_name_de.h5ad should exist
	String file_name
	
	Float? minimum_report_score
	Boolean? no_use_non_de

	Int? diskSpace = 100

	call run_scrtools_annotate {
		input:
			data_folder = data_folder,
			file_name = file_name,
			minimum_report_score = minimum_report_score,
			no_use_non_de = no_use_non_de,
			diskSpace = diskSpace
	}
}

task run_scrtools_annotate {
	String data_folder
	String file_name
	Float? minimum_report_score
	Boolean? no_use_non_de
	Int? diskSpace

	command {
		set -e
		export TMPDIR=/tmp
		python <<CODE
		from subprocess import check_call

		input_file = '${file_name}_de.h5ad'
		minimum_report_score = '${minimum_report_score}'
		no_use_non_de = '${no_use_non_de}'

		call_args = ['gsutil', '-m', 'cp', '${data_folder}/' + input_file, '.']
		print(' '.join(call_args))
		check_call(call_args)

		call_args = ['scrtools', 'annotate_cluster', input_file, '${file_name}']
		if minimum_report_score is not '':
			call_args.extend(['--minimum-report-score', minimum_report_score])
		if no_use_non_de is 'true':
			call_args.append('--do-not-use-non-de-genes')
		print(' '.join(call_args))
		check_call(call_args)
		CODE
		gsutil -m mv ${file_name}_fisher.anno.txt ${data_folder}
		gsutil -m mv ${file_name}_t.anno.txt ${data_folder}
	}

	output {
		String annotation_based_on_fisher_test = '${data_folder}/${file_name}_fisher.anno.txt'
		String annotation_based_on_t_test = '${data_folder}/${file_name}_t.anno.txt'
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
