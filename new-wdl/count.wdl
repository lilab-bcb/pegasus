workflow cellranger {
    String sampleId
    Array[String]? fastqs
    String? commaFastqs
    String referenceName
    File transcriptome_file = (if referenceName == 'GRCh38' 
                                then 'gs://fc-bcc55e6c-bec3-4b2e-9fb2-5e1526ddfcd2/reference_data/human/grch38/refdata-cellranger-GRCh38-1.2.0.tar.gz'
                                else (if referenceName == 'mm10' 
                                        then 'gs://fc-bcc55e6c-bec3-4b2e-9fb2-5e1526ddfcd2/reference_data/mouse/mm10/refdata-cellranger-mm10-1.2.0.tar.gz' 
                                        else (if referenceName == 'GRCh38_and_mm10'
                                                then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz'
                                                else referenceName)))
    Boolean? secondary
    Int? expectCells
    Int? forceCells
    Boolean? do_force_cells
    String? chemistry
    Int diskSpace
    Int memory
    Int cores
    Int preemptible

    call CellRanger {
        input:
        sampleId = sampleId,
        fastqs = fastqs,
        commaFastqs = commaFastqs,
        reference = referenceName,
        transcriptome_file = transcriptome_file,
        secondary = secondary,
        expectCells = expectCells,
        diskSpace = diskSpace,
        forceCells = forceCells,
        do_force_cells = do_force_cells,
        chemistry = chemistry,
        memory = memory,
        cores = cores,
        preemptible = preemptible
   }
   
   call ConvertCellRangerOutput {
       input:
       sampleId = sampleId,
       graphclust = CellRanger.graph_clusters,
       km2 = CellRanger.kmeans_clust_2,
       km3 = CellRanger.kmeans_clust_3,
       km4 = CellRanger.kmeans_clust_4,
       km5 = CellRanger.kmeans_clust_5,
       km6 = CellRanger.kmeans_clust_6,
       km7 = CellRanger.kmeans_clust_7,
       km8 = CellRanger.kmeans_clust_8,
       km9 = CellRanger.kmeans_clust_9,
       km10 = CellRanger.kmeans_clust_10,
       pca = CellRanger.pca_projection,
       tsne = CellRanger.tsne,
       diskSpace = diskSpace,
       memory = memory,
       cores = cores,
       preemptible = preemptible
   }

   output {
    File filtered_gene_h5 = CellRanger.filtered_gene_h5
   }
}

task CellRanger {
    String sampleId
    Array[String]? fastqs
    String? commaFastqs
    String reference
    File transcriptome_file
    Int? expectCells
    Boolean? secondary
    Int? forceCells
    Boolean? do_force_cells
    String? chemistry
    Int diskSpace
    Int memory
    Int cores
    Int preemptible

    command {
        set -e
        mkdir transcriptome_dir
        tar xf ${transcriptome_file} -C transcriptome_dir --strip-components 1
        ln -s /usr/bin/python3 python
        export PATH=$PATH:.
        python /software/scripts/orchestra_methods.py -c=count \
                                                    -id=${sampleId} \
                                                    -cf=${commaFastqs} \
                                                    -fs=["${sep='","' fastqs}"] \
                                                    -E=${expectCells} \
                                                    -F=${forceCells} \
                                                    -C=${chemistry} \
                                                    -S =${secondary}
        }
    output {
        File barcodes = "results_${sampleId}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
        File genes = "results_${sampleId}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
        File matrix = "results_${sampleId}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"
        File qc = "results_${sampleId}/outs/metrics_summary.csv"
        File report = "results_${sampleId}/outs/web_summary.html"
        File sorted_bam = "results_${sampleId}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "results_${sampleId}/outs/possorted_genome_bam.bam.bai"
        File filtered_gene_h5 = "results_${sampleId}/outs/filtered_gene_bc_matrices_h5.h5"
        File raw_gene_h5 = "results_${sampleId}/outs/raw_gene_bc_matrices_h5.h5"
        File raw_barcodes = "results_${sampleId}/outs/raw_gene_bc_matrices/${reference}/barcodes.tsv"
        File raw_genes = "results_${sampleId}/outs/raw_gene_bc_matrices/${reference}/genes.tsv"
        File raw_matrix = "results_${sampleId}/outs/raw_gene_bc_matrices/${reference}/matrix.mtx"
        File mol_info_h5 = "results_${sampleId}/outs/molecule_info.h5"
        File cloupe = "results_${sampleId}/outs/cloupe.cloupe"
        File tsne = "results_${sampleId}/outs/analysis/tsne/2_components/projection.csv"
        File graph_clusters = "results_${sampleId}/outs/analysis/clustering/graphclust/clusters.csv"
        File graph_diffexp = "results_${sampleId}/outs/analysis/diffexp/graphclust/differential_expression.csv"
        File kmeans_clust_2 = "results_${sampleId}/outs/analysis/clustering/kmeans_2_clusters/clusters.csv"
        File kmeans_clust_3 = "results_${sampleId}/outs/analysis/clustering/kmeans_3_clusters/clusters.csv" 
        File kmeans_clust_4 = "results_${sampleId}/outs/analysis/clustering/kmeans_4_clusters/clusters.csv" 
        File kmeans_clust_5 = "results_${sampleId}/outs/analysis/clustering/kmeans_5_clusters/clusters.csv" 
        File kmeans_clust_6 = "results_${sampleId}/outs/analysis/clustering/kmeans_6_clusters/clusters.csv" 
        File kmeans_clust_7 = "results_${sampleId}/outs/analysis/clustering/kmeans_7_clusters/clusters.csv" 
        File kmeans_clust_8 = "results_${sampleId}/outs/analysis/clustering/kmeans_8_clusters/clusters.csv" 
        File kmeans_clust_9 = "results_${sampleId}/outs/analysis/clustering/kmeans_9_clusters/clusters.csv" 
        File kmeans_clust_10 = "results_${sampleId}/outs/analysis/clustering/kmeans_10_clusters/clusters.csv" 
        File kmeans_diffexp_2 = "results_${sampleId}/outs/analysis/diffexp/kmeans_2_clusters/differential_expression.csv"
        File kmeans_diffexp_3 = "results_${sampleId}/outs/analysis/diffexp/kmeans_3_clusters/differential_expression.csv"
        File kmeans_diffexp_4 = "results_${sampleId}/outs/analysis/diffexp/kmeans_4_clusters/differential_expression.csv"
        File kmeans_diffexp_5 = "results_${sampleId}/outs/analysis/diffexp/kmeans_5_clusters/differential_expression.csv"
        File kmeans_diffexp_6 = "results_${sampleId}/outs/analysis/diffexp/kmeans_6_clusters/differential_expression.csv"
        File kmeans_diffexp_7 = "results_${sampleId}/outs/analysis/diffexp/kmeans_7_clusters/differential_expression.csv"
        File kmeans_diffexp_8 = "results_${sampleId}/outs/analysis/diffexp/kmeans_8_clusters/differential_expression.csv"
        File kmeans_diffexp_9 = "results_${sampleId}/outs/analysis/diffexp/kmeans_9_clusters/differential_expression.csv"
        File kmeans_diffexp_10 = "results_${sampleId}/outs/analysis/diffexp/kmeans_10_clusters/differential_expression.csv"
        File pca_components = "results_${sampleId}/outs/analysis/pca/10_components/components.csv"
        File pca_dispersion = "results_${sampleId}/outs/analysis/pca/10_components/dispersion.csv"
        File pca_genes = "results_${sampleId}/outs/analysis/pca/10_components/genes_selected.csv"
        File pca_projection = "results_${sampleId}/outs/analysis/pca/10_components/projection.csv"
        File pca_variance = "results_${sampleId}/outs/analysis/pca/10_components/variance.csv"   

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
    
task ConvertCellRangerOutput {
    String sampleId
    File graphclust
    File km2
    File km3
    File km4
    File km5
    File km6
    File km7
    File km8
    File km9
    File km10
    File pca
    File tsne
    Int diskSpace
    Int memory
    Int cores
    Int preemptible
    
    command <<<
        python /software/scripts/cell_ranger_to_scp.py \
                --graphclust ${graphclust} \
                --kmeans_2 ${km2} \
                --kmeans_3 ${km3} \
                --kmeans_4 ${km4} \
                --kmeans_5 ${km5} \
                --kmeans_6 ${km6} \
                --kmeans_7 ${km7} \
                --kmeans_8 ${km8} \
                --kmeans_9 ${km9} \
                --kmeans_10 ${km10} \
                --pca ${pca} \
                --tsne ${tsne} \
                --metadata_file "${sampleId}_metadata.txt" \
                --pca_coordinates_file "${sampleId}_pca.txt" \
                --tsne_coordinates_file "${sampleId}_tsne.txt" 
    >>>
    
    output {
        File metadata = "${sampleId}_metadata.txt"
        File pca_coords = "${sampleId}_pca.txt"
        File tsne_coords = "${sampleId}_tsne.txt"
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