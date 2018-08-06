workflow cellranger {
	# gsurl of BCL to run CellRanger mkfastq on
	String bcl
	# Csv File mapping samples, lanes and indices to bcl
	File masterCsv
	# Output directory of Fastqs
	String output_directory
	# Runtime Arguments
	Int? diskSpace = 500 
	Int? memory = 120
	Int? cores = 32
	Int? preemptible = 2

	call mkfastq {
		input:
		diskSpace =  diskSpace,
		bcl = bcl,
		output_directory = output_directory,
		masterCsv = masterCsv,
		memory = memory,
		cores = cores,
		preemptible = preemptible
	}

	output {
		String path = mkfastq.path
		String undetermined_path = mkfastq.undetermined_path
		File monitoringLog = mkfastq.monitoringLog
	}
}

task mkfastq {
	String bcl
	String output_directory
	File masterCsv
	Int diskSpace
	Int memory
	Int cores
	Int preemptible

	command {
		monitor_script.sh > monitoring.log &
		
		orchestra_methods.py -c=mkfastq \
                        -b=${bcl} \
                        -M=${masterCsv} \
                        -O=${output_directory}
	}

	output {
		String path = read_string("path.txt")
		String undetermined_path = read_string("undetermined.txt")
		File monitoringLog = "monitoring.log"
	}

	runtime {
	        docker: "singlecellportal/scrna-seq_orchestra"
	        memory: "${memory} GB"
	        bootDiskSizeGb: 12
	        disks: "local-disk ${diskSpace} HDD"
	        cpu: cores
	        preemptible: preemptible
	    }
}
