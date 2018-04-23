workflow cellranger_mkfastq_count {
    # 3 columns (Lane, Sample, Index). gs URL
	File input_csv_file
    # Sequencer output directory containing Config/, Data/, Images/, InterOp/, etc. gs URL
    String input_directory
    # Fastq output directory, gs URL
    String fastq_output_directory
    # Optional disk space for mkfastq.
    Int? mkfastq_disk_space = 500
    # 2.1.1 only currently available
    String? cell_ranger_version = "2.1.1"

	# gs URL to a transcriptome directory tar.gz
	# "gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz"
	# "gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz"
	String transcriptome

	# Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
	Boolean? secondary = false
	# If force cells, default: true
	Boolean? do_force_cells = true
	# Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. Default: 6,000 cells
	Int? force_cells = 6000
	# Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with force_cells
	Int? expect_cells = 3000
	# Optional disk space needed for cell ranger count.
	Int? count_disk_space = 250


    call cellranger_mkfastq {
        input:
            input_directory = input_directory,
            input_csv_file = input_csv_file,
            disk_space = mkfastq_disk_space,
            output_directory = fastq_output_directory
    }

	call parse_csv {
		input:
			input_csv_file = input_csv_file
	}

	scatter (sample_id in parse_csv.sample_ids) {
		call cell_ranger_count {
			input:
				sample_id = sample_id,
				data_directory = cellranger_mkfastq.output_fastqs_directory,
				transcriptome = transcriptome,
				do_force_cells = do_force_cells,
				force_cells = force_cells,
				secondary = secondary,
				expect_cells = expect_cells,
				disk_space = count_disk_space,
				cell_ranger_version = cell_ranger_version
		}
	}
}

task parse_csv {
	File input_csv_file

	command {
		set -e
		python <<CODE
		with open('${input_csv_file}') as fin:
			next(fin)
			for line in fin:
				print(line.strip().split(',')[1])
		CODE
	}

	output {
		Array[String] sample_ids = read_lines(stdout())
	}

	runtime {
		docker: "regevlab/cellranger-${cell_ranger_version}"
		preemptible: 2
	}
}

task cell_ranger_count {
	String sample_id
	String data_directory
	File transcriptome
	Boolean? do_force_cells
	Int? expect_cells
	Boolean? secondary
	Int? disk_space
	String? cell_ranger_version
	Int? force_cells

	command {
		set -e
		export TMPDIR=/tmp
		mkdir -p transcriptome_dir
		tar xf ${transcriptome} -C transcriptome_dir --strip-components 1
		gsutil -m cp -r ${data_directory}/fastqs/fastq_path/*/${sample_id} .
		python <<CODE
		from subprocess import check_call
		call_args = ['cellranger', 'count', '--id=results', '--transcriptome=transcriptome_dir', '--fastqs=${sample_id}', '--sample=${sample_id}', '--jobmode=local']
		call_args.append('--force-cells=${force_cells}' if ('${do_force_cells}' is 'true') else '--expect-cells=${expect_cells}')
		if '${secondary}' is not 'true':
			call_args.append('--nosecondary')
		check_call(call_args)
		CODE
		gsutil -m mv results/outs ${data_directory}/${sample_id}
	}

	output {
		String output_count_directory = "${data_directory}/${sample_id}"
	}

	runtime {
		docker: "regevlab/cellranger-${cell_ranger_version}"
		memory: "384 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 64
		preemptible: 2
	}
}


task cellranger_mkfastq {
	String input_directory
	File input_csv_file
	String output_directory
	Int? disk_space

	String input_name = basename(input_directory)

	command {
		set -e
		export TMPDIR=/tmp
		echo ${input_name}
		gsutil -m cp -r ${input_directory} .
		cellranger mkfastq --id=results --run=${input_name} --csv=${input_csv_file} --jobmode=local
		gsutil -m mv results/outs ${output_directory}/fastqs
	}

	output {
		String output_fastqs_directory = "${output_directory}/fastqs"
	}

	runtime {
		docker: "regevlab/cellranger-${cell_ranger_version}"
		memory: "192 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 64
		preemptible: 2
	}
}
