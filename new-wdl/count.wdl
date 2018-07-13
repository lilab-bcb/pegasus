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
    Boolean secondary = false
    Int? expectCells
    Int? forceCells
    Boolean? do_force_cells
    String? chemistry
    Int diskSpace
    Int memory
    Int cores
    Int preemptible
    File? monitoringScript = "gs://fc-6665a95b-cb65-4527-ad5e-0d1be0fdf3dc/monitor_script.sh"

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
        preemptible = preemptible,
        monitoringScript = monitoringScript
    }

    output {
        File barcodes = CellRanger.barcodes
        File genes = CellRanger.genes
        File matrix = CellRanger.matrix
        File qc = CellRanger.qc
        File report = CellRanger.report
        File sorted_bam = CellRanger.sorted_bam
        File sorted_bam_index = CellRanger.sorted_bam_index
        File filtered_gene_h5 = CellRanger.filtered_gene_h5
        File raw_gene_h5 = CellRanger.raw_gene_h5
        File raw_barcodes = CellRanger.raw_barcodes
        File raw_genes = CellRanger.raw_genes
        File raw_matrix = CellRanger.raw_matrix
        File mol_info_h5 = CellRanger.mol_info_h5
        File cloupe = CellRanger.cloupe
        File monitoringLog = CellRanger.monitoringLog
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
    File monitoringScript

    command {
        chmod u+x ${monitoringScript}
        ${monitoringScript} > monitoring.log &

        orchestra_methods.py -c=count \
                            -id=${sampleId} \
                            -cf=${commaFastqs} \
                            -fs "${sep='" "' fastqs}" \
                            -E=${expectCells} \
                            -F=${forceCells} \
                            -C=${chemistry} \
                            -S=${secondary} \
                            -tf=${transcriptome_file}
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
        File monitoringLog = "monitoring.log"
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