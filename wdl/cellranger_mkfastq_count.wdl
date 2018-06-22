workflow cellranger_mkfastq_count {
	# 4 columns (Sample, Flowcell, Lane, Index). gs URL
	File input_csv_file
	# CellRanger output directory, gs URL
	String cellranger_output_directory
	# Optional disk space for mkfastq.
	Int? mkfastq_disk_space = 1000
	# 2.1.1 only currently available
	String? cellranger_version = "2.1.1"
	# Whether to delete input_directory. If false, you should delete this folder yourself so as to not incur storage charges.
	Boolean? delete_input_directory = false

	# gs://regevlab-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz
	# gs://regevlab-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz
	# gs://regevlab-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz
	# GRCh38, mm10, GRCh38_and_mm10, or a URL to a tar.gz file
	String transcriptome
	
	File transcriptome_file = (if transcriptome == 'GRCh38' 
								then 'gs://regevlab-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz'
								else (if transcriptome == 'mm10' 
										then 'gs://regevlab-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz' 
										else (if transcriptome == 'GRCh38_and_mm10'
												then 'gs://regevlab-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz'
												else transcriptome)))

	# Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
	Boolean? secondary = false
	# If force cells, default: true
	Boolean? do_force_cells = true
	# Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. Default: 6,000 cells
	Int? force_cells = 6000
	# Expected number of recovered cells. Default: 3,000 cells. Mutually exclusive with force_cells
	Int? expect_cells = 3000
	# Optional disk space needed for cell ranger count.
	Int? count_disk_space = 500

	# Number of cpus per cellranger job
	Int? num_cpu = 64
	# Number of preemptible tries 
	Int? preemptible = 2

	call parse_input_csv {
		input:
			input_csv_file = input_csv_file,
			output_dir = cellranger_output_directory,
			cellranger_version = cellranger_version,
			preemptible = preemptible
	}

	scatter (run_id in parse_input_csv.run_ids) {
		call cellranger_mkfastq {
			input:
				run_id = run_id,
				input_directory = parse_input_csv.inpdirs[run_id],
				input_csv_file = parse_input_csv.bcls[run_id],
				output_directory = cellranger_output_directory,
				disk_space = mkfastq_disk_space,
				num_cpu = num_cpu,
				delete_input_directory = delete_input_directory,
				cellranger_version = cellranger_version,
				preemptible = preemptible
		}
	}

	scatter (sample_id in parse_input_csv.sample_ids) {
		call cellranger_count {
			input:
				sample_id = sample_id,
				run_ids = parse_input_csv.sample2run[sample_id],
				data_directory = cellranger_output_directory,
				transcriptome = transcriptome_file,
				do_force_cells = do_force_cells,
				force_cells = force_cells,
				secondary = secondary,
				expect_cells = expect_cells,
				disk_space = count_disk_space,
				num_cpu = num_cpu,
				cellranger_version = cellranger_version,
				preemptible = preemptible,
				input_fastqs_directory = cellranger_mkfastq.output_fastqs_directory
		}
	}

	call collect_summaries {
		input:
			summaries = cellranger_count.output_metrics_summary,
			sample_ids = cellranger_count.output_count_directory,
			cellranger_version = cellranger_version,
			preemptible = preemptible
	}
}

task parse_input_csv {
	File input_csv_file
	String output_dir
	String? cellranger_version
	Int? preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE

		import os
		import pandas as pd 
		from subprocess import check_call
		from collections import defaultdict

		df = pd.read_csv('${input_csv_file}', header = 0)
		sample2run = defaultdict(list)
		with open('run_ids.txt', 'w') as fo1, open('inpdirs.txt', 'w') as fo2, open('bcls.txt', 'w') as fo3:
			for input_dir in df['Flowcell'].unique():
				run_id = os.path.basename(input_dir)
				bcl_df = df.loc[df['Flowcell'] == input_dir, ['Lane','Sample', 'Index']]
				bcl_df.to_csv(run_id + '_bcl.csv', index = False)
				for sample_id in bcl_df['Sample']:
					sample2run[sample_id].append(run_id)

				# copy bcl CSV to google bucket
				call_args = ['gsutil', '-q', 'cp', run_id + '_bcl.csv', '${output_dir}/']
				check_call(call_args)

				fo1.write(run_id + '\n')
				fo2.write(run_id + '\t' + input_dir + '\n')
				fo3.write(run_id + '\t${output_dir}/' + run_id + '_bcl.csv\n')

		with open('sample_ids.txt', 'w') as fo1, open('sample2run.txt', 'w') as fo2:
			for key, value in sample2run.items():
				fo1.write(key + '\n')
				fo2.write(key + '\t' + ','.join(value) + '\n')

		CODE
	}

	output {
		Array[String] run_ids = read_lines('run_ids.txt')
		Map[String, String] inpdirs = read_map('inpdirs.txt')
		Map[String, String] bcls = read_map('bcls.txt')
		Array[String] sample_ids = read_lines('sample_ids.txt')
		Map[String, String] sample2run = read_map('sample2run.txt')
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		preemptible: "${preemptible}"
	}
}

task cellranger_mkfastq {
	String run_id
	String input_directory
	File input_csv_file
	String output_directory
	Int? disk_space
	Int? num_cpu
	String? cellranger_version
	Boolean? delete_input_directory
	Int? preemptible

	String input_name = basename(input_directory)

	command {
		set -e
		export TMPDIR=/tmp
		gsutil -q -m cp -r ${input_directory} .
		cellranger mkfastq --id=results --run=${input_name} --csv=${input_csv_file} --jobmode=local
		gsutil -q -m mv results/outs ${output_directory}/${run_id}
		
		python <<CODE

		from subprocess import check_call, check_output, CalledProcessError

		if '${delete_input_directory}' is 'true':
			try:
				call_args = ['gsutil', '-q', 'stat', '${output_directory}/${run_id}/qc_summary.json']
				check_output(call_args)
				call_args = ['gsutil', '-q' 'rm' '-r', '${input_directory}']
				check_call(call_args)
			except CalledProcessError:
				print("Failed to move outputs to Google bucket.")
		
		CODE
	}

	output {
		String output_fastqs_directory = "${output_directory}/${run_id}"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		memory: "192 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}

task cellranger_count {
	String sample_id
	String run_ids
	String data_directory
	File transcriptome
	Boolean? do_force_cells
	Int? expect_cells
	Boolean? secondary
	Int? disk_space
	Int? num_cpu
	String? cellranger_version
	Int? force_cells
	Int? preemptible
	Array[String] input_fastqs_directory

	command {
		set -e
		export TMPDIR=/tmp
		mkdir -p transcriptome_dir
		tar xf ${transcriptome} -C transcriptome_dir --strip-components 1

		python <<CODE

		from subprocess import check_call

		fastqs = []
		for run_id in '${run_ids}'.split(','):
			call_args = ['gsutil', '-q', '-m', 'cp', '-r', '${data_directory}' + '/' + run_id + '/fastq_path/*/' + '${sample_id}', '.']
			check_call(call_args)
			call_args = ['mv', '${sample_id}', '${sample_id}_' + run_id]
			check_call(call_args)
			fastqs.append('${sample_id}_' + run_id)
	
		call_args = ['cellranger', 'count', '--id=results', '--transcriptome=transcriptome_dir', '--fastqs=' + ','.join(fastqs), '--sample=${sample_id}', '--jobmode=local']
		call_args.append('--force-cells=${force_cells}' if ('${do_force_cells}' is 'true') else '--expect-cells=${expect_cells}')
		if '${secondary}' is not 'true':
			call_args.append('--nosecondary')
		check_call(call_args)

		CODE

		gsutil -q -m cp -r results/outs ${data_directory}/${sample_id}
	}

	output {
		String output_count_directory = "${data_directory}/${sample_id}"
		File output_metrics_summary = "results/outs/metrics_summary.csv"
		File output_web_summary = "results/outs/web_summary.html"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		memory: "384 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
	}
}

task collect_summaries {
	Array[File] summaries
	Array[String] sample_ids
	String? cellranger_version
	Int? preemptible

	command {
		python <<CODE

		import pandas as pd
		import os
		import xlsxwriter

		summaries = pd.read_csv('${write_lines(summaries)}', header = None)
		sample_ids = pd.read_csv('${write_lines(sample_ids)}', header = None).applymap(lambda x: os.path.basename(x))
		df_info = pd.concat([summaries, sample_ids], axis = 1)
		df_info.columns = ['summary', 'sample_id']
		dfs = []
		for idx, row in df_info.iterrows():
			df = pd.read_csv(row['summary'], header = 0)
			df.index = [row['sample_id']]
			dfs.append(df)
		df_res = pd.concat(dfs)
		df_res.index.name = "Sample"
		writer = pd.ExcelWriter('summaries.xlsx', engine = 'xlsxwriter')
		df_res.to_excel(writer, sheet_name = "summaries")
		writer.save()

		CODE
	}

	output {
		File metrics_summaries = "summaries.xlsx"
	}

	runtime {
		docker: "regevlab/cellranger-${cellranger_version}"
		preemptible: "${preemptible}"
	}
}
