workflow cellranger {
	String diskSpace
	File bcl
	File masterCsv

	call getCSV as csvTask{
		input:
		diskSpace =  diskSpace,
		masterCsv = masterCsv
	}

	call getFlowcell as fl{
		input:
		bcl = bcl,
		diskSpace = diskSpace,
		csv = csvTask.csv,
	}

	call CellRanger{
		input:
		bcl = bcl,
		diskSpace = diskSpace,
		csv = csvTask.csv,
		flowcell = fl.flowcell
	}

	output{
		Array[File] sampleFastQs = CellRanger.sampleFastQs
		File qc = CellRanger.qc
	}
}

task getCSV{
	String diskSpace
	File masterCsv

	command {
		# create the CSV File for this flowcell
		set -e
		mkdir fastqs
		ln -s /usr/bin/python3 python
		export PATH=$PATH:.
		python <<CODE
		import os
		import pandas as pd
		master_csv = "${masterCsv}"
		df = pd.read_csv(master_csv,header=0)
		#at this point leaving no sorting, for when we move to cloud this will work
		#df = df.loc[df['flowcell'] == "mkfastq"]
		df = df.drop(columns=['flowcell'])
		df.to_csv('sample.csv',index=False)
		CODE
	}
	output {
		File csv = "sample.csv"
	}

	runtime {
	    docker: "singlecellportal/cell-ranger-count-2.1.1"
	    memory: "16 GB"
	    bootDiskSizeGb: 12
	    disks: "local-disk ${diskSpace} HDD"
	    cpu: 1
	    preemptible: 2
	}
}

task getFlowcell {
	String diskSpace
	File bcl
	File csv

		command {
	set -e
	export PATH=$PATH:.
	python <<CODE
	import os
	import tarfile
	import xml.etree.ElementTree as ET
	from subprocess import call
	bcl = '${bcl}'
	csv = '${csv}'
	try:
		tarfile.is_tarfile(bcl)
		tar = tarfile.open(bcl)
		tar.extractall()
		names = tar.getnames()
		tar.close()
		run = names[0].replace("./._", "./")
	except:
		print("Not a Tar Error")

	# get flowcell

	tree = ET.parse(run+"/RunInfo.xml")
	root = tree.getroot()
	flowcell = root.find('Run').find('Flowcell').text
	print(flowcell)
	CODE

	}

		output{
			String flowcell = read_string(stdout())
		}
		runtime {
		    docker: "singlecellportal/cell-ranger-count-2.1.1"
		    memory: "16 GB"
		    bootDiskSizeGb: 12
		    disks: "local-disk ${diskSpace} HDD"
		    cpu: 1
		    preemptible: 2
		}

}

task CellRanger{
	String diskSpace
	File bcl
	File csv
	String flowcell

	command {
set -e
mkdir fastqs
ln -s /usr/bin/python3 python
export PATH=$PATH:.
python <<CODE
import os
import tarfile
import xml.etree.ElementTree as ET
from subprocess import call
bcl = '${bcl}'
csv = '${csv}'
try:
	tarfile.is_tarfile(bcl)
	tar = tarfile.open(bcl)
	tar.extractall()
	names = tar.getnames()
	tar.close()
	run = names[0].replace("./._", "./")
except:
	print("Not a Tar Error")

# get flowcell

tree = ET.parse(run+"/RunInfo.xml")
root = tree.getroot()
flowcell = root.find('Run').find('Flowcell').text
print(flowcell)

call_args = list()
call_args.append('cellranger')
call_args.append('mkfastq')
call_args.append('--run=' + run)
call_args.append('--csv=' + csv)
call_args.append('--output-dir=fastqs')
call(call_args)

# export and tar fastqs
CODE

cd fastqs/${flowcell}
python <<CODE	
import os
import tarfile
from subprocess import call
samples = os.listdir()
for sample in samples:
	call_args = list()
	call_args.append('tar')
	call_args.append('cf')
	call_args.append(sample + '.tar')
	call_args.append(sample)
	call(call_args)
CODE
}

	output{
		Array[File] sampleFastQs = glob("fastqs/H35KCBCXY/*.tar")
		File qc = "H35KCBCXY/outs/qc_summary.json"
	}
	runtime {
	    docker: "singlecellportal/cell-ranger-count-2.1.1"
	    memory: "16 GB"
	    bootDiskSizeGb: 12
	    disks: "local-disk ${diskSpace} HDD"
	    cpu: 1
	    preemptible: 2
	}
}