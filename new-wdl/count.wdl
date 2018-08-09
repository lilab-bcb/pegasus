workflow cellranger {
    # Name of Sample to run CellRanger Count on
    String sampleId
    # Array of Fastq Directories for running decoupled
    Array[String]? fastqs
    # Comma Seperated String of Fastq Directories for Running in orchestration
    String? commaFastqs
    # Reference Name for Count
    String reference
    # Reference File for Count
    File transcriptome_file
    # Run CellRanger built in secondary analysis? Default to false because we have Bo's pipeline
    Boolean secondary = false
    # CellRanger count input
    Int? expectCells
    # CellRanger count input
    Int? forceCells
    # CellRanger count input
    Boolean? do_force_cells
    # CellRanger count input
    String? chemistry
    # Runtime Arguments
    Int? diskSpace = 500
    Int? memory = 120
    Int? cores = 32
    Int? preemptible = 2

     call CellRanger {
        input:
        sampleId = sampleId,
        fastqs = fastqs,
        commaFastqs = commaFastqs,
        transcriptome_file = transcriptome_file,
        reference = reference,
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

    command {
        monitor_script.sh > monitoring.log &

        orchestra_methods.py -c=count \
                            -id=${sampleId} \
                            -cf=${commaFastqs} \
                            -fs "${sep='" "' fastqs}" \
                            -E=${expectCells} \
                            -F=${forceCells} \
                            -C=${chemistry} \
                            -S=${secondary} \
                            -tf=${transcriptome_file} \
                            -dfc=${do_force_cells}
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
        File monitoringLog = "monitoring.log"
    }
    
 runtime {
         docker: "singlecellportal/scrna-seq_orchestra"
         memory: "${memory} GB"
         bootDiskSizeGb: 12
         disks: "local-disk ${diskSpace} HDD"
         cpu: cores
         preemptible: preemptible
     }
}