workflow cellranger_mkfastq {
	# Input BCL directory, gs url
	String input_bcl_directory
	# 3 column CSV file (Lane, Sample, Index)
	File input_csv_file
	# CellRanger output directory, gs url
	String output_directory

	# Whether to delete input bcl directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean delete_input_bcl_directory
	# Currently, only 2.1.1 is available
	String? cellranger_version = "2.1.1"

	# Number of cpus per cellranger job
	Int? num_cpu = 64
	# Memory in GB
	Int? memory = 128
	# Disk space in GB
	Int? disk_space = 1500
	# Number of preemptible tries 
	Int? preemptible = 2

	call run_cellranger_mkfastq {
		input:
			input_bcl_directory = input_bcl_directory,
			input_csv_file = input_csv_file,
			output_directory = output_directory,
			delete_input_bcl_directory = delete_input_bcl_directory,
			cellranger_version = cellranger_version,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible
	}

	output {
		String output_fastqs_directory = run_cellranger_mkfastq.output_fastqs_directory
		String output_fastqs_flowcell_directory = run_cellranger_mkfastq.output_fastqs_flowcell_directory
		File monitoringLog = run_cellranger_mkfastq.monitoringLog
	}
}

task run_cellranger_mkfastq {
	String input_bcl_directory
	File input_csv_file
	String output_directory
	Boolean delete_input_bcl_directory
	String cellranger_version
	Int num_cpu
	Int memory
	Int disk_space
	Int preemptible

	String run_id = basename(input_bcl_directory)

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &
		gsutil -q -m cp -r ${input_bcl_directory} .
		# cp -r ${input_bcl_directory} .
		cellranger mkfastq --id=results --run=${run_id} --csv=${input_csv_file} --jobmode=local
		gsutil -q -m rsync -d -r results/outs ${output_directory}/${run_id}_fastqs
		# cp -r results/outs ${output_directory}/${run_id}_fastqs


		python <<CODE
		import os
		from subprocess import check_call, check_output, CalledProcessError
		with open("output_fastqs_flowcell_directory.txt", "w") as fout:
			flowcell = [name for name in os.listdir('results/outs/fastq_path') if name != 'Reports' and name != 'Stats' and os.path.isdir('results/outs/fastq_path/' + name)][0]
			fout.write('${output_directory}/${run_id}_fastqs/fastq_path/' + flowcell + '\n')
		if '${delete_input_bcl_directory}' is 'true':
			try:
				call_args = ['gsutil', '-q', 'stat', '${output_directory}/${run_id}_fastqs/qc_summary.json']
				print(' '.join(call_args))
				check_output(call_args)
				call_args = ['gsutil', '-q', '-m', 'rm', '-r', '${input_bcl_directory}']
				print(' '.join(call_args))
				check_call(call_args)
				print('${input_bcl_directory} is deleted!')
			except CalledProcessError:
				print("Failed to move outputs to Google bucket.")
		CODE
	}

	output {
		String output_fastqs_directory = "${output_directory}/${run_id}_fastqs"
		String output_fastqs_flowcell_directory = select_first(read_lines("output_fastqs_flowcell_directory.txt"))
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		memory: "${memory} GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}
