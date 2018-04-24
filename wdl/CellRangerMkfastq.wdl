workflow cellranger_mkfastq {
	String input_directory
	File input_csv_file
	String output_directory
	Int? diskSpace = 500
	Int? numCPU = 64

	call run_cellranger_mkfastq {
		input:
			input_directory = input_directory,
			input_csv_file = input_csv_file,
			output_directory = output_directory,
			diskSpace = diskSpace,
			numCPU = numCPU
	}
}

task run_cellranger_mkfastq {
	String input_directory
	File input_csv_file
	String output_directory
	Int? diskSpace
	Int? numCPU

	String input_name = basename(input_directory)

	command {
		set -e
		export TMPDIR=/tmp
		gsutil -m cp -r ${input_directory} .
		cellranger mkfastq --id=results --run=${input_name} --csv=${input_csv_file} --jobmode=local
		gsutil -m mv results/outs ${output_directory}/fastqs
	}

	output {
		String output_fastqs_directory = "${output_directory}/fastqs"
	}

	runtime {
		docker: "regevlab/cellranger_2.1.1"
		memory: "192 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: "${numCPU}"
		preemptible: 2
	}
}
