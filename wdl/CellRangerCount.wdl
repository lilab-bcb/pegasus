workflow cellranger_count {
	File input_csv_file
	# A gs link, assume fastqs is under data_directory; We will also move all cellranger outputs into the data_directory
	String data_directory

	# gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz
	# gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz
	# mm10, GRCh38, or a URL to a tar.gz file
	File transcriptome

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

	call parse_csv {
		input:
			input_csv_file = input_csv_file
	}
    if(transcriptome=='mm10') {
        transcriptome = 'gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz'
    }
    if(transcriptome=='GRCh38') {
        transcriptome = 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz'
    }
	scatter (sampleId in parse_csv.sampleIds) {
		call CellRangerCount {
			input:
				sampleId = sampleId,
				data_directory = data_directory,
				transcriptome = transcriptome,
				do_force_cells = do_force_cells,
				forceCells = forceCells,
				secondary = secondary,
				expectCells = expectCells,
				diskSpace = diskSpace,
				version = version
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
		Array[String] sampleIds = read_lines(stdout())
	}

	runtime {
		docker: "regevlab/cellranger-${version}"
		preemptible: 2
	}
}

task CellRangerCount {
	String sampleId
	String data_directory
	File transcriptome
	Boolean? do_force_cells 
	Int? expectCells
	Boolean? secondary
	Int? diskSpace
	String? version
	Int? forceCells

	command {
		set -e
		export TMPDIR=/tmp
		mkdir -p transcriptome_dir
		tar xf ${transcriptome} -C transcriptome_dir --strip-components 1
		gsutil -m cp -r ${data_directory}/fastqs/fastq_path/*/${sampleId} .
		python <<CODE
		from subprocess import check_call
		call_args = ['cellranger', 'count', '--id=results', '--transcriptome=transcriptome_dir', '--fastqs=${sampleId}', '--sample=${sampleId}', '--jobmode=local']
		call_args.append('--force-cells=${forceCells}' if ('${do_force_cells}' is 'true') else '--expect-cells=${expectCells}')
		if '${secondary}' is not 'true':
			call_args.append('--nosecondary')
		check_call(call_args)
		CODE
		gsutil -m mv results/outs ${data_directory}/${sampleId}
	}

	output {
		String output_count_directory = "${data_directory}/${sampleId}"
	}

	runtime {
		docker: "regevlab/cellranger-${version}"
		memory: "384 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${diskSpace} HDD"
		cpu: 64
		preemptible: 2
	}
}
