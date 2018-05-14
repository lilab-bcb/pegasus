workflow cellranger_mkfastq_count {
	# 3 columns (Lane, Sample, Index). gs URL
	File input_csv_file
	# Sequencer output directory containing Config/, Data/, Images/, InterOp/, etc. gs URL
	String input_directory
	# CellRanger output directory, gs URL
	String cellranger_output_directory
	# Optional disk space for mkfastq.
	Int? mkfastq_disk_space = 500
	# 2.1.1 only currently available
	String? cellranger_version = "2.1.1"
	# Whether to delete input_directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean? delete_input_directory = false

	# gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz
	# gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz
	# mm10, GRCh38, or a URL to a tar.gz file
	String transcriptome
	
	File transcriptome_file = (if transcriptome == 'GRCh38' 
								then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz'
								else (if transcriptome == 'mm10' 
										then 'gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz' 
										else transcriptome))

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

	# Number of cpus per cellranger job
	Int? num_cpu = 64


	call cellranger_mkfastq {
		input:
			input_directory = input_directory,
			input_csv_file = input_csv_file,
			output_directory = cellranger_output_directory,
			disk_space = mkfastq_disk_space,
			num_cpu = num_cpu,
			delete_input_directory=delete_input_directory,
			cellranger_version =cellranger_version
	}

	scatter (sample_id in cellranger_mkfastq.sample_ids) {
		call cellranger_count {
			input:
				sample_id = sample_id,
				data_directory = cellranger_output_directory,
				transcriptome = transcriptome_file,
				do_force_cells = do_force_cells,
				force_cells = force_cells,
				secondary = secondary,
				expect_cells = expect_cells,
				disk_space = count_disk_space,
				num_cpu = num_cpu,
				cellranger_version = cellranger_version
		}
	}
}

task cellranger_mkfastq {
	String input_directory
	File input_csv_file
	String output_directory
	Int? disk_space
	Int? num_cpu
	String? cellranger_version
	Boolean? delete_input_directory

	String input_name = basename(input_directory)

	command {
		set -e
		export TMPDIR=/tmp
		python <<CODE
		import os
		with open('${input_csv_file}') as fin, open('sample_ids.txt', 'w') as fout:
			next(fin)
			for line in fin:
				fout.write(line.strip().split(',')[1] + '\n')
		CODE
		gsutil -q -m cp -r ${input_directory} .
		cellranger mkfastq --id=results --run=${input_name} --csv=${input_csv_file} --jobmode=local
		gsutil -q -m mv results/outs ${output_directory}/fastqs
		if delete_input_directory and os.path.exists("${output_directory}/fastqs"):
		    gsutil -q rm -r ${input_directory}
	}

	output {
		String output_fastqs_directory = "${output_directory}/fastqs"
		Array[String] sample_ids = read_lines("sample_ids.txt")
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		memory: "192 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: 2
	}
}

task cellranger_count {
	String sample_id
	String data_directory
	File transcriptome
	Boolean? do_force_cells
	Int? expect_cells
	Boolean? secondary
	Int? disk_space
	Int? num_cpu
	String? cellranger_version
	Int? force_cells

	command {
		set -e
		export TMPDIR=/tmp
		mkdir -p transcriptome_dir
		tar xf ${transcriptome} -C transcriptome_dir --strip-components 1
		gsutil -q -m cp -r ${data_directory}/fastqs/fastq_path/*/${sample_id} .
		python <<CODE
		from subprocess import check_call
		call_args = ['cellranger', 'count', '--id=results', '--transcriptome=transcriptome_dir', '--fastqs=${sample_id}', '--sample=${sample_id}', '--jobmode=local']
		call_args.append('--force-cells=${force_cells}' if ('${do_force_cells}' is 'true') else '--expect-cells=${expect_cells}')
		if '${secondary}' is not 'true':
			call_args.append('--nosecondary')
		check_call(call_args)
		CODE
		gsutil -q -m mv results/outs ${data_directory}/${sample_id}
	}

	output {
		String output_count_directory = "${data_directory}/${sample_id}"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		memory: "384 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: 2
	}
}
