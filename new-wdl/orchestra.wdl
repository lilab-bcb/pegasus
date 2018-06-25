import "scRNA-Seq/new-wdl/mkfastq.wdl" as mkfastq
import "scRNA-Seq/new-wdl/count.wdl" as count

workflow orchestra {
  String id
  File run
  File csv
  String flowcell
  String sampleId
  String referenceName
  File transcriptomeTarGz
  Boolean? secondary
  Int? expectCells
  String diskSpace
  Int? forceCells
  Boolean? do_force_cells
  String? chemistry
	
	call mkfastq.CellRanger as fq {
		input:
		id = id,
		run = run,
		csv = csv,
		diskSpace = diskSpace,
		flowcell = flowcell
	}

	call count.CellRanger as cnt {
		input:
   		sampleId = sampleId,
   		fastqs = fq.sampleFastQs,
   		reference = referenceName,
   		transcriptomeTarGz = transcriptomeTarGz,
   		secondary = secondary,
   		expectCells = expectCells,
   		diskSpace = diskSpace,
   		forceCells = forceCells,
   		do_force_cells = do_force_cells,
   		chemistry = chemistry
    }

    call count.ConvertCellRangerOutput {
    	input:
    	sampleId = sampleId,
    	graphclust = cnt.graph_clusters,
    	km2 = cnt.kmeans_clust_2,
    	km3 = cnt.kmeans_clust_3,
    	km4 = cnt.kmeans_clust_4,
    	km5 = cnt.kmeans_clust_5,
    	km6 = cnt.kmeans_clust_6,
    	km7 = cnt.kmeans_clust_7,
    	km8 = cnt.kmeans_clust_8,
    	km9 = cnt.kmeans_clust_9,
    	km10 = cnt.kmeans_clust_10,
    	pca = cnt.pca_projection,
    	tsne = cnt.tsne,
    	diskSpace = diskSpace
    }
}