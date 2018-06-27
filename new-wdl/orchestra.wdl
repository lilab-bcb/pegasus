import "scRNA-Seq/new-wdl/mkfastq.wdl" as mkfastq
import "scRNA-Seq/new-wdl/count.wdl" as count

workflow orchestra {
  Array[File] bcls
  Array[String] sampleIds
  File masterCsv
  String referenceName
  File transcriptomeTarGz
  Boolean? secondary
  Int? expectCells
  String diskSpace
  Int? forceCells
  Boolean? do_force_cells
  String? chemistry
	
  # scatter per flow cell
  scatter(bcl in bcls){
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
  scatter(sampleId in sampleIds){
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

