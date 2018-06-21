workflow cellranger_count {
	File input_csv_file
	# A gs link. We assume each sample has a subfolder containing FASTQs under input_directory
	String input_directory
	# A gs link. Each sample will have a subfolder containing 'cellranger count' outputs under output_directory
	String output_directory

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
	# Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expectCells. Default: 6,000 cells
	Int? forceCells = 6000
	# Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with forceCells
	Int? expectCells = 3000
	# 2.1.1 or 2.0.2
	String? version = "2.1.1"
	# Disk space needed per task. If you don't know enter "250"
	Int? diskSpace = 250
	# Number of cpus per cellranger count job
	Int? numCPU = 64

	call parse_csv {
		input:
			input_csv_file = input_csv_file,
			version = version
	}

	scatter (sampleId in parse_csv.sampleIds) {
		call CellRangerCount {
			input:
				sampleId = sampleId,
				input_directory = input_directory,
				output_directory = output_directory,
				transcriptome = transcriptome_file,
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
	String? version
	
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
		Array[String] sampleIds = read_lines(stdout())
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

	command {
		set -e
		export TMPDIR=/tmp
		mkdir -p transcriptome_dir
		tar xf ${transcriptome} -C transcriptome_dir --strip-components 1
		gsutil -m cp -r ${input_directory}/${sampleId} .
		python <<CODE
		from subprocess import check_call
		call_args = ['cellranger', 'count', '--id=results', '--transcriptome=transcriptome_dir', '--fastqs=${sampleId}', '--sample=${sampleId}', '--jobmode=local']
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
