workflow cellranger {
	String diskSpace
	File run
	File csv
	String flowcell
	call CellRanger{
		input:
		id = id,
		run = run,
		csv = csv,
		diskSpace = diskSpace,
		flowcell = flowcell
	}

	output{
		Array[File] sampleFastQs = CellRanger.sampleFastQs
		File sampleCSV = csv
		File qc = CellRanger.qc
	}
}

task CellRanger{
	String diskSpace
	File run
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
from subprocess import call
run = '${run}'
csv = '${csv}'
flowcell = '${flowcell}'
try:
	tarfile.is_tarfile(run)
	tar = tarfile.open(run)
	tar.extractall()
	names = tar.getnames()
	tar.close()
	run = names[0].replace("./._", "./")
except:
	print(run, " is not a Tar")
call_args = list()
call_args.append('cellranger')
call_args.append('mkfastq')
call_args.append('--run=' + run)
call_args.append('--csv=' + csv)
call_args.append('--output-dir=fastqs')
call(call_args)


CODE
cd fastqs/${flowcell}
python <<CODE
import os
import tarfile
from subprocess import call
# Compress every sample output folder
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
		Array[File] sampleFastQs = glob("fastqs/${flowcell}/*.tar")
		File qc = "${flowcell}/outs/qc_summary.json"
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