import "scRNA-Seq/new-wdl/mkfastq.wdl" as mkfastq
import "scRNA-Seq/new-wdl/count.wdl" as count

workflow orchestra {
  File masterCsv
  String referenceName
  File transcriptomeTarGz
  Boolean? secondary
  Int? expectCells
  String diskSpace
  Int? forceCells
  Boolean? do_force_cells
  String? chemistry
	
  call parseCsv as parse {
    input:
      masterCsv = masterCsv
      diskSpace = diskSpace
  }

  # scatter per flow cell
  scatter(bcl in parse.bcls){
      call mkfastq.cellranger as fq {
        input:
          bcl = bcl,
          masterCsv = masterCsv,
          diskSpace = diskSpace
      }
  }
	
  # gathered- a bunch of fast qs for every sample across flowcells
  # scatter by sample name
  # in count 
  scatter(sampleId in parse.sampleIds){
      call count.cellranger as cnt {
        input:
          sampleId = sampleId,
          fastqs = flatten(fq.sampleFastQs),
          referenceName = referenceName,
          transcriptomeTarGz = transcriptomeTarGz,
          secondary = secondary,
          expectCells = expectCells,
          diskSpace = diskSpace,
          forceCells = forceCells,
          do_force_cells = do_force_cells,
          chemistry = chemistry
      }
    }	
}

task parseCsv {
  File masterCsv
  String diskSpace

  command {
    set -e
    mkdir fastqs
    ln -s /usr/bin/python3 python
    export PATH=$PATH:.
    python <<CODE
    import os
    import pandas as pd
    master_csv = "${masterCsv}"
    df = pd.read_csv(master_csv,header=0)
    bcl_paths = set(df['flowcell'].tolist())
    bcl_file = open('bcls.txt', 'w+')
    for item in bcl_paths:
      bcl_file.write("%s\n" % item)
    bcl_file.close()
    sampleIds = set(df['sample'].tolist())
    samples_file = open('samples.txt', 'w+')
    for item in sampleIds:
      samples_file.write("%s\n" % item)
    samples_file.close()
    CODE
  }

  output {
    Array[String] bcls = read_lines('bcls.txt')
    Array[String] sampleIds = read_lines('samples.txt')
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

