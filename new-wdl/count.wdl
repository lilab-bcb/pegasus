workflow cellranger {
    String sampleId
    Array[String] fastqs
    String referenceName
    File transcriptome_file = (if referenceName == 'GRCh38' 
                                then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz'
                                else (if referenceName == 'mm10' 
                                        then 'gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz' 
                                        else (if referenceName == 'GRCh38_and_mm10'
                                                then 'gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38_and_mm10-1.2.0.tar.gz'
                                                else referenceName)))
    Boolean? secondary
    Int? expectCells
    String diskSpace
    Int? forceCells
    Boolean? do_force_cells
    String? chemistry
    String memory
    String cores
    String preemptible
    String countOutPath

    call CellRanger {
        input:
        sampleId = sampleId,
        fastqs = fastqs,
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
        countOutPath = countOutPath
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
       preemptible = preemptible,
       countOutPath = countOutPath
   }
}

task CellRanger {
    String sampleId
    Array[String] fastqs
    String reference
    File transcriptome_file
    Int? expectCells
    Boolean? secondary
    String diskSpace
    Int? forceCells
    Boolean? do_force_cells
    String? chemistry
    String memory
    String cores
    String preemptible
    String countOutPath

    command {
        set -e
        mkdir transcriptome_dir
        tar xf ${transcriptome_file} -C transcriptome_dir --strip-components 1
        ln -s /usr/bin/python3 python
        export PATH=$PATH:.
        python <<CODE
        import os
        import tarfile
        from subprocess import call
        dirs = set()
        sample = "${sampleId}"
        i = 0
        for f in ["${sep='","' fastqs}"]:
            # get the fastqs
            call(["mkdir", str(i)])
            call(["gsutil", "-q", "-m", "cp", "-r", f, str(i)])
            dirs.add(str(i)+"/" +sample+"/")
            i+=1
        expect_cells = '${expectCells}'
        force_cells = '${forceCells}'
        secondary = '${secondary}'
        chemistry = '${chemistry}'
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--sample=' + sample)
        call_args.append('--id=results_'+sample)
        call_args.append('--fastqs=' + ','.join(dirs))
        if secondary is not 'true':
            call_args.append('--nosecondary')
        if force_cells is not '':
            call_args.append('--force-cells=' + str(force_cells))
        elif expect_cells is not '':
            call_args.append('--expect-cells=' + str(expect_cells))
        if chemistry is not '':
            call_args.append('--chemistry='+chemistry)
        call(call_args)

        # Move the outputs to the output directory
        call(["gsutil", "-q", "-m", "cp", "-r", "results_"+sample+"/", "${countOutPath}"])
        CODE
        }
    output {
        String out_path = "${countOutPath}results_${sampleId}"
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
    String memory
    String cores
    String preemptible
    String countOutPath
    
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
        gsutil -q -m cp -r "${sampleId}_metadata.txt" "${countOutPath}results_${sampleId}"
        gsutil -q -m cp -r "${sampleId}_pca.txt" "${countOutPath}results_${sampleId}" 
        gsutil -q -m cp -r "${sampleId}_tsne.txt" "${countOutPath}results_${sampleId}"  
    >>>
    
    output {
        String out_path = "${countOutPath}results_${sampleId}"
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