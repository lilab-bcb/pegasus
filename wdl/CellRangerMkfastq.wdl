workflow cellranger_mkfastq {
	File input_directory_tar
	File input_csv_file
	String output_directory
	Int? diskSpace = 500

	call run_cellranger_mkfastq {
		input:
			input_directory_tar = input_directory_tar,
			input_csv_file = input_csv_file,
			diskSpace = diskSpace,
			output_directory = output_directory
	}
}

task run_cellranger_mkfastq {
	File input_directory_tar
	File input_csv_file
	String output_directory
	Int diskSpace

	command {
		export TMPDIR=/tmp
		mkdir -p illumina_data
		tar -xf ${input_directory_tar} -C illumina_data
		cellranger mkfastq --id=results --run=illumina_data --csv=${input_csv_file} --jobmode=local
		gsutil -m mv results/outs ${output_directory}/fastqs
		echo "Done"
	}

	output {
		String status = stdout()
	}

	runtime {
		docker: "bigbadbo/cellranger_2.1.1"
		memory: "192 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: 64
		preemptible: 2
	}
}
