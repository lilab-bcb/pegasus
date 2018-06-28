workflow cellranger {
	String diskSpace
	String bcl
	File masterCsv
	String output_directory
	String memory
	String cores
	String preemptible

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
	String diskSpace
	String output_directory
	File masterCsv
	String memory
	String cores
	String preemptible

	command {
		set -e
		export PATH=$PATH:.
		gsutil -q -m cp -r -L log.txt ${bcl} .
		mkdir fastqs
		ln -s /usr/bin/python3 python
		python <<CODE
		import os
		import tarfile
		import xml.etree.ElementTree as ET
		import pandas as pd
		from subprocess import call
		# open up log.txt
		df = pd.read_csv('log.txt', header=0)
		path = "${bcl}"
		run = list(filter(None, path.split("/")))[-1]

		# get flowcell
		tree = ET.parse(run+"/RunInfo.xml")
		root = tree.getroot()
		flowcell = root.find('Run').find('Flowcell').text

		# Create the csv
		master_csv = "${masterCsv}"
		df = pd.read_csv(master_csv,header=0)
		df = df.loc[df['Flowcell'] == path]
		df = df.drop(columns=['Flowcell'])
		df.to_csv('sample.csv',index=False)

		call_args = list()
		call_args.append('cellranger')
		call_args.append('mkfastq')
		call_args.append('--run=' + run+"/")
		call_args.append('--csv=' + 'sample.csv')
		call_args.append('--output-dir=fastqs')
		call(call_args)

		# output the fastqs
		call(["gsutil", "-q", "-m", "mv", 'fastqs/'+flowcell, "${output_directory}"])
		call(["gsutil", "-q", "-m", "mv", flowcell+"/outs/qc_summary.json", "${output_directory}"])
		
		file = open("path.txt","w") 
		file.write("${output_directory}"+flowcell+"/")
		file.close()
		CODE
	}

	output {
		String path = read_string("path.txt")
	}

	runtime {
	    docker: "singlecellportal/cell-ranger-count-2.1.1"
	    memory: "${memory} GB"
	    bootDiskSizeGb: 12
	    disks: "local-disk ${diskSpace} HDD"
	    cpu: "${cores}"
	    preemptible: "${preemptible}"
	}
}
