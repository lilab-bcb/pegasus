workflow cellranger {
    String sampleId
    Array[File] fastqs
    String referenceName
    File transcriptomeTarGz
    Boolean? secondary
    Int? expectCells
    String diskSpace
    Int? forceCells
    Boolean? do_force_cells

    call CellRanger {
        input:
        sampleId = sampleId,
        fastqs = fastqs,
        reference = referenceName,
        transcriptomeTarGz = transcriptomeTarGz,
        secondary = secondary,
        expectCells = expectCells,
        diskSpace = diskSpace,
        forceCells = forceCells,
        do_force_cells = do_force_cells
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
       diskSpace = diskSpace
   }
}

task CellRanger {
    String sampleId
    Array[File] fastqs
    String reference
    File transcriptomeTarGz
    Int? expectCells
    Boolean? secondary
    String diskSpace
    Int? forceCells
    Boolean? do_force_cells

    command {
        set -e
        mkdir transcriptome_dir
        tar xf ${transcriptomeTarGz} -C transcriptome_dir --strip-components 1
        ln -s /usr/bin/python3 python
        export PATH=$PATH:.
        python <<CODE
        import os
        import tarfile
        from subprocess import call
        dirs = dict()
        for f in ["${sep='","' fastqs}"]:
            file = os.path.dirname(f)
            try:
                tarfile.is_tarfile(f)
                tar = tarfile.open(f)
                tar.extractall()
                names = tar.getnames()
                tar.close()
                dirs.setdefault(names[0], True)
            except:
                dirs.setdefault(file, True)
        expect_cells = '${expectCells}'
        force_cells = '${forceCells}'
        secondary = '${secondary}'
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--id=' + '${sampleId}')
        call_args.append('--fastqs=' + ','.join(list(dirs.keys())))
        if secondary is not 'true':
            call_args.append('--nosecondary')
        if force_cells is not '':
            call_args.append('--force-cells=' + str(force_cells))
        elif expect_cells is not '':
            call_args.append('--expect-cells=' + str(expect_cells))
        call(call_args)
        CODE
        }
    output {
        File barcodes = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
        File genes = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
        File matrix = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"
        File qc = "${sampleId}/outs/metrics_summary.csv"
        File report = "${sampleId}/outs/web_summary.html"
        File sorted_bam = "${sampleId}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "${sampleId}/outs/possorted_genome_bam.bam.bai"
        File filtered_gene_h5 = "${sampleId}/outs/filtered_gene_bc_matrices_h5.h5"
        File raw_gene_h5 = "${sampleId}/outs/raw_gene_bc_matrices_h5.h5"
        File raw_barcodes = "${sampleId}/outs/raw_gene_bc_matrices/${reference}/barcodes.tsv"
        File raw_genes = "${sampleId}/outs/raw_gene_bc_matrices/${reference}/genes.tsv"
        File raw_matrix = "${sampleId}/outs/raw_gene_bc_matrices/${reference}/matrix.mtx"
        File mol_info_h5 = "${sampleId}/outs/molecule_info.h5"
        File cloupe = "${sampleId}/outs/cloupe.cloupe"
        File tsne = "${sampleId}/outs/analysis/tsne/2_components/projection.csv"
        File graph_clusters = "${sampleId}/outs/analysis/clustering/graphclust/clusters.csv"
        File graph_diffexp = "${sampleId}/outs/analysis/diffexp/graphclust/differential_expression.csv"
        File kmeans_clust_2 = "${sampleId}/outs/analysis/clustering/kmeans_2_clusters/clusters.csv"
        File kmeans_clust_3 = "${sampleId}/outs/analysis/clustering/kmeans_3_clusters/clusters.csv" 
        File kmeans_clust_4 = "${sampleId}/outs/analysis/clustering/kmeans_4_clusters/clusters.csv" 
        File kmeans_clust_5 = "${sampleId}/outs/analysis/clustering/kmeans_5_clusters/clusters.csv" 
        File kmeans_clust_6 = "${sampleId}/outs/analysis/clustering/kmeans_6_clusters/clusters.csv" 
        File kmeans_clust_7 = "${sampleId}/outs/analysis/clustering/kmeans_7_clusters/clusters.csv" 
        File kmeans_clust_8 = "${sampleId}/outs/analysis/clustering/kmeans_8_clusters/clusters.csv" 
        File kmeans_clust_9 = "${sampleId}/outs/analysis/clustering/kmeans_9_clusters/clusters.csv" 
        File kmeans_clust_10 = "${sampleId}/outs/analysis/clustering/kmeans_10_clusters/clusters.csv" 
        File kmeans_diffexp_2 = "${sampleId}/outs/analysis/diffexp/kmeans_2_clusters/differential_expression.csv"
        File kmeans_diffexp_3 = "${sampleId}/outs/analysis/diffexp/kmeans_3_clusters/differential_expression.csv"
        File kmeans_diffexp_4 = "${sampleId}/outs/analysis/diffexp/kmeans_4_clusters/differential_expression.csv"
        File kmeans_diffexp_5 = "${sampleId}/outs/analysis/diffexp/kmeans_5_clusters/differential_expression.csv"
        File kmeans_diffexp_6 = "${sampleId}/outs/analysis/diffexp/kmeans_6_clusters/differential_expression.csv"
        File kmeans_diffexp_7 = "${sampleId}/outs/analysis/diffexp/kmeans_7_clusters/differential_expression.csv"
        File kmeans_diffexp_8 = "${sampleId}/outs/analysis/diffexp/kmeans_8_clusters/differential_expression.csv"
        File kmeans_diffexp_9 = "${sampleId}/outs/analysis/diffexp/kmeans_9_clusters/differential_expression.csv"
        File kmeans_diffexp_10 = "${sampleId}/outs/analysis/diffexp/kmeans_10_clusters/differential_expression.csv"
        File pca_components = "${sampleId}/outs/analysis/pca/10_components/components.csv"
        File pca_dispersion = "${sampleId}/outs/analysis/pca/10_components/dispersion.csv"
        File pca_genes = "${sampleId}/outs/analysis/pca/10_components/genes_selected.csv"
        File pca_projection = "${sampleId}/outs/analysis/pca/10_components/projection.csv"
        File pca_variance = "${sampleId}/outs/analysis/pca/10_components/variance.csv"   
    }
    
    runtime {
        docker: "singlecellportal/cell-ranger-count-2.0.2"
        memory: "416 GB"
        bootDiskSizeGb: 12
        disks: "local-disk ${diskSpace} HDD"
        cpu: 64
        preemptible: 2
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
    String diskSpace
    
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
                --tsne_coordinates_file "${sampleId}_tsne.txt" \
    >>>
    
    output {
        File metadata = "${sampleId}_metadata.txt"
        File pca_coords = "${sampleId}_pca.txt"
        File tsne_coords = "${sampleId}_tsne.txt"
    }

    runtime {
        docker: "singlecellportal/cell-ranger-count-2.0.2"
        memory: "416 GB"
        bootDiskSizeGb: 12
        disks: "local-disk ${diskSpace} HDD"
        cpu: 64
        preemptible: 2
    }
}