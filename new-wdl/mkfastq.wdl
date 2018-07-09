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
		set -e
		export PATH=$PATH:.
		gsutil -q -m cp -r -L log.txt ${bcl} .
		mkdir fastqs
		ln -s /usr/bin/python3 python
		python /software/scripts/orchestra_methods.py -c=mkfastq \
                                                    -b=${bcl} \
                                                    -M=${masterCsv} \
                                                    -O=${output_directory}
	}

	output {
		String path = read_string("path.txt")
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
