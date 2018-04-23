# scRNA-Seq

Cell Ranger mkfastq/count inputs:

Name | Description | Example | Default
--- | --- | --- | ---
input_csv_file | # 3 columns (Lane, Sample, Index) | "gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/my_file.csv"
input_directory | Sequencer output directory containing Config/, Data/, Images/, InterOp/, etc. | "gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/my_dir" | 
fastq_output_directory | Fastq output directory | "gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/my_dir" |
mkfastq_disk_space | Optional disk space for mkfastq. | 500 | 500 
version | Cell ranger version | "2.1.1" | "2.1.1"
transcriptome | mm10, GRCh38, or a gs URL to a transcriptome directory tar.gz | "mm10" | 
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
do_force_cells | force cells | true | true
force_cells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expect_cells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
count_disk_space | Disk space needed for cell ranger count | 500 | 250

