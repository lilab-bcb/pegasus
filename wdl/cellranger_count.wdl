workflow cellranger_count {
	File input_csv_file
	# A gs link. We assume each sample has a subfolder containing FASTQs under input_directory
	String input_directory
	# A gs link. Each sample will have a subfolder containing 'cellranger count' outputs under output_directory
	String output_directory

	# gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz
	# gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz
	# gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz
	# GRCh38, mm10, GRCh38_and_mm10, or a URL to a tar.gz file
	String transcriptome
	
	File transcriptome_file = (if transcriptome == 'GRCh38' 
								then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz'
								else (if transcriptome == 'mm10' 
										then 'gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz' 
										else (if transcriptome == 'GRCh38_and_mm10'
												then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz'
												else transcriptome)))

	# Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
	Boolean? secondary = false
	# If force cells, default: true
	Boolean? do_force_cells = true
	# Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expectCells. Default: 6,000 cells
	Int? forceCells = 6000
	# Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with forceCells
	Int? expectCells = 3000
	# 2.1.1 or 2.0.2
	String? version = "2.1.1"
	# Disk space needed per task. If you don't know enter "250"
	Int? diskSpace = 500
	# Number of cpus per cellranger count job
	Int? numCPU = 64
	# Comma-separated list of chemistries for each sample
	String? chemistry_string

	call parse_csv {
		input:
			input_csv_file = input_csv_file,
			chemistry_str = chemistry_string,
			version = version
	}

	scatter (i in range(length(parse_csv.sampleIds))) {
		call CellRangerCount {
			input:
				sampleId = parse_csv.sampleIds[i],
				input_directory = input_directory,
				output_directory = output_directory,
				transcriptome = transcriptome_file,
				chemistry = parse_csv.chemistry[i],
				do_force_cells = do_force_cells,
				forceCells = forceCells,
				secondary = secondary,
				expectCells = expectCells,
				diskSpace = diskSpace,
				numCPU = numCPU,
				version = version
		}
	}
}

task parse_csv {
	File input_csv_file
	String? chemistry_str
	String? version
	
	command {
		set -e
		python <<CODE
		n = 0
		with open('${input_csv_file}') as fin:
			next(fin)
			for line in fin:
				print(line.strip().split(',')[1])
				n += 1

		with open('chemistry.txt', 'w') as fout:
			my_str = '${chemistry_str}'
			if my_str == '':
				my_str = ','.join(['auto'] * n)
			fout.write('\n'.join(my_str.split(',')) + '\n')
		CODE
	}

	output {
		Array[String] sampleIds = read_lines(stdout())
		Array[String] chemistry = read_lines('chemistry.txt')
	}

	runtime {
		docker: "regevlab/cellranger-${version}"
		preemptible: 2
	}
}

task CellRangerCount {
	String sampleId
	String input_directory
	String output_directory
	File transcriptome
	Boolean? do_force_cells 
	Int? expectCells
	Boolean? secondary
	Int? diskSpace
	Int? numCPU
	String? version
	Int? forceCells
	String? chemistry = "auto"

	command {
		set -e
		export TMPDIR=/tmp
		mkdir -p transcriptome_dir
		tar xf ${transcriptome} -C transcriptome_dir --strip-components 1
		gsutil -m cp -r ${input_directory}/${sampleId} .
		python <<CODE
		from subprocess import check_call
		call_args = ['cellranger', 'count', '--id=results', '--transcriptome=transcriptome_dir', '--fastqs=${sampleId}', '--sample=${sampleId}', '--jobmode=local', '--chemistry=${chemistry}']
		call_args.append('--force-cells=${forceCells}' if ('${do_force_cells}' is 'true') else '--expect-cells=${expectCells}')
		if '${secondary}' is not 'true':
			call_args.append('--nosecondary')
		check_call(call_args)
		CODE
		gsutil -m mv results/outs ${output_directory}/${sampleId}
	}

	output {
		String output_count_directory = "${output_directory}/${sampleId}"
	}

	runtime {
		docker: "regevlab/cellranger-${version}"
		memory: "384 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: "${numCPU}"
		preemptible: 2
	}
}
