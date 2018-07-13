workflow cellranger {
	String bcl
	File masterCsv
	String output_directory
	Int diskSpace
	Int memory
	Int cores
	Int preemptible

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
		orchestra_methods.py -c=mkfastq \
                        -b=${bcl} \
                        -M=${masterCsv} \
                        -O=${output_directory}
	}

	output {
		String path = read_string("path.txt")
		String undetermined_path = read_string("undetermined.txt")
	}

	runtime {
	        docker: "singlecellportal/cell-ranger-count-2.1.1"
	        memory: "${memory} GB"
	        bootDiskSizeGb: 12
	        disks: "local-disk ${diskSpace} HDD"
	        cpu: cores
	        preemptible: preemptible
	    }
}
