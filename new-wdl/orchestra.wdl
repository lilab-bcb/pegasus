import "scRNA-Seq/new-wdl/mkfastq.wdl" as mkfastq
import "scRNA-Seq/new-wdl/count.wdl" as count

workflow orchestra {
  File masterCsv
  String referenceName
  String fastq_output_directory
  File transcriptomeTarGz
  Boolean? secondary
  Int? expectCells
  String diskSpace
  Int? forceCells
  Boolean? do_force_cells
  String? chemistry
  String memory
  String cores
  String preemptible
	
  call parseCsv as parse {
    input:
      masterCsv = masterCsv,
      diskSpace = diskSpace,
      preemptible = preemptible,
      cores = cores,
      memory = memory
  }

  # scatter per flow cell
  scatter(bcl in parse.bcls){
      call mkfastq.cellranger as fq {
        input:
          bcl = bcl,
          masterCsv = masterCsv,
          diskSpace = diskSpace,
          output_directory = fastq_output_directory,
          preemptible = preemptible,
          cores = cores,
          memory = memory
      }
  }
	
  # gathered- a bunch of fast qs for every sample across flowcells
  # scatter by sample name
  # in count 
  scatter(sampleId in parse.sampleIds){
      call filterSamples as filter {
        input:
          paths = fq.path,
          sample = sampleId,
          diskSpace = diskSpace,
          preemptible = preemptible,
          cores = cores,
          memory = memory
      }

      call count.cellranger as cnt {
        input:
          sampleId = sampleId,
          fastqs = filter.filteredPaths,
          referenceName = referenceName,
          transcriptomeTarGz = transcriptomeTarGz,
          secondary = secondary,
          expectCells = expectCells,
          diskSpace = diskSpace,
          forceCells = forceCells,
          do_force_cells = do_force_cells,
          chemistry = chemistry,
          preemptible = preemptible,
          cores = cores,
          memory = memory
      }
    }	
}

task parseCsv {
  File masterCsv
  String diskSpace
  String memory
  String cores
  String preemptible

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
    bcl_paths = set(df['Flowcell'].tolist())
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
        memory: "${memory} GB"
        bootDiskSizeGb: 12
        disks: "local-disk ${diskSpace} HDD"
        cpu: "${cores}"
        preemptible: "${preemptible}"
    }

}

task filterSamples {
  Array[String] paths
  String sample
  String diskSpace
  String memory
  String cores
  String preemptible

  command {
    set -e
    mkdir fastqs
    ln -s /usr/bin/python3 python
    export PATH=$PATH:.
    python <<CODE
    import os
    import pandas as pd
    from subproccess import run

    # Now we have a list of every sample
    samples = []
    for f in ["${sep='","' paths}"]:
      samples = samples + run(["gsutil", "ls", f]).split("\n")

    key = "${sample}"
    paths_list = open('paths.txt', 'w+')
    for path in samples:
      if path.endswith(key):
        bcl_file.write("%s\n" % path)
    paths_list.close()
  }

  output {
    Array[String] filteredPaths = read_lines("paths.txt")
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

