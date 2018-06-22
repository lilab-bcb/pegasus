# scRNA-Seq

Table of Contents
-----------------

* [First Time Running - Authenticate with Google](#first_time)
* [Run Cell Ranger mkfastq/count](#run_cellranger)
  * [Cell Ranger mkfastq/count inputs](#cellranger_input)
  * [Cell Ranger mkfastq/count outputs](#cellranger_output)
* [Run Cell Ranger count only](#run_cellranger_count)
  * [Cell Ranger count inputs](#cellranger_count_input)
* [Run Single Cell RNA-Seq analysis tools](#run_scrtools)
  * [Prepare count_matrix.csv](#count_matrix_csv)
  * [scrtools aggregate_matrix](#aggr_mat)
  * [scrtools merge_matrix](#merge_mat)
  * [scrtools cluster](#cluster)
  * [scrtools annotate](#annotate)
  * [scrtools subcluster](#subcluster)

## <a name="first_time"></a> First Time Running - Authenticate with Google

1. Ensure the Google Cloud SDK is installed. 
    
    Note: Broad users do not have to install this-they can type *use Google-Cloud-SDK* to make the Google Cloud tools available. 

1. Type *gcloud auth login* to login to Google Cloud.

1. Copy and paste the link in your unix terminal into your web browser.

1. Enter authorization code in unix terminal


## <a name="run_cellranger"></a> Run Cell Ranger mkfastq/count
1. [Create](https://software.broadinstitute.org/firecloud/documentation/article?id=10746) a [FireCloud](https://portal.firecloud.org/) workspace or use an existing workspace.

1. Copy your sequencing output to the workspace using gsutil in your unix terminal. 

    It is highly recommended that you delete the **BCL files** after the pipeline is finished by turning on the **delete_input_directory** option.

    Example: *gsutil -m rsync -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4*
    
    -m means copy in parallel, -r means copy the directory recursively.
    
    Note: Broad users need to be on an UGER node (not a login node) in order to use the -m flag
    
    You can also read [FireCloud instructions](https://software.broadinstitute.org/firecloud/documentation/article?id=10574) on uploading data.

1. Create a sample sheet. A sample sheet is a simple CSV with sample, flowcell, lane, and index columns, which describe the way to demultiplex the flowcells. The flowcell column should be google bucket urls obtained in the previous step. The index column should contain a 10x sample set name (e.g., SI-GA-A12). If a sample is sequenced in multiple flowcells, this smaple should be listed multiple times in the CSV (one line per flowcell). In the following example, we have 4 samples sequenced in two flowcells.
   
   Example: 
    ```
    Sample,Flowcell,Lane,Index
    sample_1,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8
    sample_2,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8
    sample_3,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8
    sample_4,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8
    sample_1,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8
    sample_2,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8
    sample_3,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8
    sample_4,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8
    ```
1. Upload your sample sheet to the workspace.  

    Example: *gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

1. Import cellranger_mkfastq_count method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type cellranger_mkfastq_count.

1. Uncheck "Configure inputs/outputs using the Workspace Data Model"

### <a name="cellranger_input"></a> Cell Ranger mkfastq/count inputs:

*Cell Ranger mkfastq/count* takes Illumina outputs as input and runs *cellranger mkfastq* and *cellranger count*. Please see the description of inputs below. Note that required inputs are shown in bold.

Name | Description | Example | Default
--- | --- | --- | ---
**input\_csv\_file** | Sample Sheet (contains Sample, Flowcell, Lane, Index) | "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv" | 
**transcriptome** | keyword (**GRCh38** for human, **mm10** for mouse, **GRCh38\_and\_mm10** for human and mouse) or gs URL to a transcriptome tar.gz | "GRCh38" | 
**cellranger\_output\_directory** | Cellranger output directory | "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output" |
**delete\_input\_directory** | If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges. | true | 
cellranger_version | Cellranger version | "2.1.1" | "2.1.1"
do_force_cells | force cells | true | true
force_cells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expect_cells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
mkfastq_disk_space | Optional disk space in gigabytes for mkfastq. | 1000 | 1000 
count_disk_space | Disk space in gigabytes needed for cell ranger count | 500 | 500
num_cpu | Number of cpus to request for one node | 64 | 64
preemptible | Number of preemptible tries | 2 | 2

### <a name="cellranger_output"></a> Cell Ranger mkfastq/count outputs:

See the table below for important *Cell Ranger mkfastq/count* outputs.

Name | Type | Description
--- | --- | ---
output_fastqs_directory | Array[String] | A list of google bucket urls containing FASTQ files, one url per flowcell.
output_count_directory | Array[String] | A list of google bucket urls containing count matrices, one url per sample.
metrics_summaries | File | A excel spreadsheet containing QCs for each sample.
output_web_summary | Array[File] | A list of htmls visualizing QCs for each sample (cellranger count output).

## <a name="run_cellranger_count"></a> Run Cell Ranger count only

Sometimes, people might want to perform demux locally and run **cellranger count** on the cloud. This section describes how to perform counting on cloud.

1. [Create](https://software.broadinstitute.org/firecloud/documentation/article?id=10746) a [FireCloud](https://portal.firecloud.org/) workspace or use an existing workspace.

1. Copy your FASTQ files to the workspace using gsutil in your unix terminal. 

    You should upload a folder, which contains one subfolder per sample. Each subfolder contains all FASTQ files from the corresponding sample.

    Example: *gsutil -m rsync -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4*
    
    -m means copy in parallel, -r means copy the directory recursively.
    
    Note: Broad users need to be on an UGER node (not a login node) in order to use the -m flag
    
    You can also read [FireCloud instructions](https://software.broadinstitute.org/firecloud/documentation/article?id=10574) on uploading data.

1. Create a sample sheet. A sample sheet is a simple CSV with lane, sample, and index columns. You should have this CSV file for running **cellranger mkfastq** locally. 
   
   Example: 
    ```
    Lane,Sample,Index
    *,sample_1,SI-GA-A8
    *,sample_2,SI-GA-B8
    *,sample_3,SI-GA-C8
    *,sample_4,SI-GA-D8
    ```
1. Upload your sample sheet to the workspace.  

    Example: *gsutil cp /foo/bar/projects/my_bcl.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

1. Import cellranger_count method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type cellranger_count.

1. Uncheck "Configure inputs/outputs using the Workspace Data Model"

### <a name="cellranger_count_input"></a> Cell Ranger count inputs:

Name | Description | Example | Default
--- | --- | --- | ---
**input\_csv\_file** | Sample Sheet (contains Lane, Sample, Index) | "gs://fc-e0000000-0000-0000-0000-000000000000/my_file.csv" | 
**input\_directory** | gs URL containing FASTQ files | "gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4" | 
**transcriptome** | keyword (**GRCh38** for human, **mm10** for mouse, **GRCh38\_and\_mm10** for human and mouse) or gs URL to a transcriptome tar.gz | "GRCh38" | 
**output\_directory** | Cellranger count output directory | "gs://fc-e0000000-0000-0000-0000-000000000000/my_dir" |
do_force_cells | force cells | true | true
forceCells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expectCells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
numCPU | Number of cpus to request for one node | 64 | 64
diskSpace | Disk space in gigabytes needed for cell ranger count | 500 | 500
version | Cellranger version | "2.1.1" | "2.1.1"

### <a name="cellranger_count_output"></a> Cell Ranger count outputs:

See the table below for important *Cell Ranger count* outputs.

Name | Type | Description
--- | --- | ---
output_count_directory | Array[String] | A list of google bucket urls containing count matrices, one url per sample.

## <a name="run_scrtools"></a> Run Single Cell RNA-Seq analysis tools (scrtools)

Before we run the scrtools, we need to first prepare a CSV file, *count_matrix.csv*, which describes the metadata for each 10x channel.

### <a name="count_matrix_csv"></a> count_matrix.csv

```
Sample,Source,Platform,Donor,Reference,Location

S1,bone_marrow,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
S2,bone_marrow,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
S3,pbmc,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
S4,pbmc,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
```

Please note that *Sample*, *Reference*, and *Location* are required. *Sample* refers to the sample names listed in the BCL CSV file, *Reference* refers to the genome name, and *Location* refers to the cellranger output folder. The cellranger count output file filtered_gene_bc_matrices_h5.h5 is expected to be located at *Location/Sample* (e.g. gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S1). You are free to add any other columns and these columns will be used in selecting channels for futher analysis. In this example, we have *Source*, which refers to the tissue of origin, *Platform*, which refers to the sequencing platform, and *Donor*, which refers to donor ID.

You should upload **count_matrix.csv** to your workspace

    Example: *gsutil cp /foo/bar/projects/my_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

Then you can aggregate 10x count matrices into a single count matrix using **scrtools aggregate_matrix**

### <a name="aggr_mat"></a> scrtools aggregate_matrix

*scrtools aggregate_matrix* is used to aggregate individual 10x channels into a big count matrix for downstream analysis. Please see the inputs below.

Name | Description | Example | Default
--- | --- | --- | ---
**input_count_matrix_csv** | Input CSV file describing metadata of each 10x channel | "my_count_matrix.csv" | 
**output_folder** | This is the folder for all analysis results | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**output_name** | Output file name of this task, the count matrix *output_name_10x.h5* and metadata file *output_name.attr.csv* will be generated | "my_aggr_mat_bm" | 
genome | The genome cellranger used to generate count matrices | "GRCh38" | "GRCh38"
restrictions | Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. Multiple restrictions are separated by ';'. If not restrictions are provided, all channels in the *count_matrix.csv* will be selected | "Source:bone_marrow" | 
attributes | Specify a comma-separated list of outputted attributes (i.e. column names in *count_matrix.csv*). These attributes will be imported to *output_name.attr.csv* | "Source,Donor" | 
groupby | When we know there are different groups in the study, such as bone_marrow and pbmc, we could perform batch correction separately for each group. This optional field is used to create a new attribute, GroupBy, that helps to identify groups. It takes the format of 'value' or 'attribute=value' or 'attr1+attr2+....' | "Source+Donor" | 
diskSpace | Disk space in gigabytes needed for this task | 100 | 100

### <a name="merge_mat"></a> scrtools merge_matrix

*scrtools merge_matrix* is used to combine two or more aggregated matrices produced by *scrtools aggregate_matrix* into one big count matrix. This step is optional. 

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**input_names** | A comma-separated list of input names for count matrices that you want to merge | "my_aggr_mat_bm,my_aggr_mat_pbmc" | 
**output_name** | Output file name of this task, the count matrix *output_name_10x.h5* and metadata file *output_name.attr.csv* will be generated | "my_aggr_mat" | 
genome | The genome cellranger used to generate count matrices | "GRCh38" | "GRCh38"
symbols | A comma-separated list of symbols representing each input matrix, this input is used with 'attributes' | "bm,pbmc" | 
attributes | A comma-separated list of attributes. When merging matrices, the matrix symbol defined in 'symbols' will be added in front of these attributes | "Donor" | 
diskSpace | Disk space in gigabytes needed for this task | 100 | 100

### <a name="cluster"></a> scrtools cluster

*scrtools cluster* performs PCA, tSNE visualization, diffusion map, Louvain clustering, and differential expression analysis (fisher's exact test and t test). This is the main step people want to know. Please see inputs below.


Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**input_name** | Input name of the aggreagated/merged matric | "my_aggr_mat" | 
**output_name** | Output file name of this task | "results" | 
num_cpu | Number of CPUs to use | 64 | 64
genome | The genome cellranger used to generate count matrices | "GRCh38" | "GRCh38"
output_filtration_results | A boolean value indicating if you want to output QC metrics into a spreadsheet | true | false
correct_batch_effect | A boolean value indicating if you want to correct for batch effects | true | false
batch_group_by | If correct for batch effects, if you want to correct it separately for each group. Group is defined by this input. It should be an attribute (e.g. GroupBy) | "Source" | 
plot_by_side | By default, this step will generate one tSNE plot colored by louvain cluster labels. If you provide an attribute name here, it will put the same tSNE colored by the attribute on the right-hand side | "Source" | 
legend_on_data | A boolean value indicating if we should put legends on data | true | false
plot_composition | A boolean value indicating if you want to generate composition plots for the attribute in *plot_by_side*, which for each cluster shows the percentage of cells in each attribute value | true | false
figure_size | Sizes of the composition plots in inches | "10,8" | "6,4"
plot_diffusion_map | A boolean value indicating if interactive 3D diffusion maps should be generated | true | false 
de_analysis | A boolean value indicating if you want to perform differential expression analysis | true | true
fold_change | Minimum fold change in either percentage (fisher test) or log expression (t test) to report a DE gene in spreadsheets | 2 | 1.5
labels | Cluster labels that help de_analysis identify clusters | "louvain_labels" | "louvain_labels"
output_loom | A boolean value indicating if you want to output loom-formatted files | false | false
import_attributes | Import attributes contained in the comma-separated list into the analysis object | "Donor" | 
min_genes | The minimum number of expressed genes to be consider a valid cell | 500 | 500
max_genes | *max_genes* - 1 is the maximum number of expressed genes to be consider a valid cell | 6000 | 6000
mito_prefix | Prefix of mitochondrial genes | "MT-" | "MT-"
percent_mito | Only keep cells with mitochondrial ratio less than *percent_mito* | 0.1 | 0.1
gene_percent_cells | Genes expressed in with less than *gene_percent_cells* * *number_of_cells* cells will be excluded from variable gene selection step. The default value requires a gene expressed in at least 3 cells when there are 6,000 cells | 0.0005 | 0.0005
counts_per_cell_after | Normalize each cell so that the sum of its normalized counts is equal to *counts_per_cell_after* | 1e5 | 1e5
louvain_resolution | Resolution parameter of the louvain clustering algorithm | 1 | 1.3
diskSpace | Disk space in gigabytes needed for this task | 250 | 250

Here are the description of outputs.

Name | Description | Required output
--- | --- | --- 
output_name.tsne.png | tSNE plots colored by louvain cluster labels and user-specified attribute | Yes
output_name.h5ad | Results in Scanpy's anndata format | Yes
output_name_var.h5ad | Results containing only variable genes in Scanpy's anndata format | Yes
output_name.filt.xlsx | Spreadsheet containing number of cells kept for each channel after filtering | No, only present if *output_filtration_results* is set
output_name.composition.frequency.png | Composition plot. Each cluster has a stacked bar showing the frequency of cells from each condition in this cluster | No, only present if *plot_composition* is set
output_name.composition.normalized.png | Computation plot. Each cluster has non-stacked bars showing the percentage of cells within each condition that belong to this cluster | No, only present if *plot_composition* is set
output_name.diffmap_cluster.html | Interactive 3D diffusion map colored by louvain cluster labels | No, only present if *plot_diffusion_map* is set
output_name.diffmap_condition.html | Interactive 3D diffusion map colored by the attribute in *plot_by_side* | No, only present if both *plot_diffusion_map* and *plot_by_side* are set
output_name_de.h5ad | Differential expression results stored in Scanpy's anndata format | No, only present if *de_analysis* is set
output_name_de_analysis_fisher.xlsx | Spreadsheet shows up-regulated and down-regulated genes for each cluster, calculated by Fisher's exact test | No, only present if *de_analysis* is set
output_name_de_analysis_t.xlsx | Spreadsheet shows up-regulated and down-regulated genes for each cluster, calculated by T test | No, only present if *de_analysis* is set
output_name.loom | Results in loom format | No, only present if *output_loom* is set
output_name_var.loom | Results containing only variable genes in loom format | No, only present if *output_loom* is set

### <a name="annotate"></a> scrtools annotate

This step outputs putative cell type annotations for each cluster based on known markers. It required differential expression analyses performed in *scrtools cluster* step. This step is optional and currently it only works for human immune cells. Please see the inputs below.

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**file_name** | Input file name, *file_name_de.h5ad* should exist | "results" | 
minimum_report_score | This step calculates a score between [0, 1] for each pair of cluster and putative cell type. It only report putative cell types with a minimum score of *minimum_report_score* | 0.1 | 0.5
no_use_non_de | Do not consider non-differentially-expressed genes as down-regulated | true | false
diskSpace | Disk space in gigabytes needed for this task | 100 | 100

### <a name="subcluster"></a> scrtools subcluster

This step performs subcluster analysis based on *scrtools cluster* outputs. Please see the inputs below.

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**input_name** | Input name of *scrtools cluster* result | "results" | 
**output_name** | Output file name of this task | "results_sub" | 
**cluster_ids** | A comma-separated list of cluster IDs (numbered from 1) for subcluster | "1,3,9" | 
num_cpu | Number of CPUs to use | 64 | 64
correct_batch_effect | A boolean value indicating if you want to correct for batch effects | true | false
plot_by_side | By default, this step will generate one tSNE plot colored by louvain cluster labels. If you provide an attribute name here, it will put the same tSNE colored by the attribute on the right-hand side | "Source" | 
legend_on_data | A boolean value indicating if we should put legends on data | true | false
plot_composition | A boolean value indicating if you want to generate composition plots for the attribute in *plot_by_side*, which for each cluster shows the percentage of cells in each attribute value | true | false
figure_size | Sizes of the composition plots in inches | "10,8" | "6,4"
plot_diffusion_map | A boolean value indicating if interactive 3D diffusion maps should be generated | true | false 
de_analysis | A boolean value indicating if you want to perform differential expression analysis | true | true
fold_change | Minimum fold change in either percentage (fisher test) or log expression (t test) to report a DE gene in spreadsheets | 2 | 1.5
labels | Cluster labels that help de_analysis identify clusters | "louvain_labels" | "louvain_labels"
louvain_resolution | Resolution parameter of the louvain clustering algorithm | 1 | 1.3
output_loom | A boolean value indicating if you want to output loom-formatted files | false | false
diskSpace | Disk space in gigabytes needed for this task | 250 | 250

Here are the outputs.

Name | Description | Required output
--- | --- | --- 
output_name.tsne.png | tSNE plots colored by louvain cluster labels and user-specified attribute | Yes
output_name.h5ad | Results in Scanpy's anndata format | Yes
output_name_var.h5ad | Results containing only variable genes in Scanpy's anndata format | Yes
output_name.filt.xlsx | Spreadsheet containing number of cells kept for each channel after filtering | No, only present if *output_filtration_results* is set
output_name.composition.frequency.png | Composition plot. Each cluster has a stacked bar showing the frequency of cells from each condition in this cluster | No, only present if *plot_composition* is set
output_name.composition.normalized.png | Computation plot. Each cluster has non-stacked bars showing the percentage of cells within each condition that belong to this cluster | No, only present if *plot_composition* is set
output_name.diffmap_cluster.html | Interactive 3D diffusion map colored by louvain cluster labels | No, only present if *plot_diffusion_map* is set
output_name.diffmap_condition.html | Interactive 3D diffusion map colored by the attribute in *plot_by_side* | No, only present if both *plot_diffusion_map* and *plot_by_side* are set
output_name_de.h5ad | Differential expression results stored in Scanpy's anndata format | No, only present if *de_analysis* is set
output_name_de_analysis_fisher.xlsx | Spreadsheet shows up-regulated and down-regulated genes for each cluster, calculated by Fisher's exact test | No, only present if *de_analysis* is set
output_name_de_analysis_t.xlsx | Spreadsheet shows up-regulated and down-regulated genes for each cluster, calculated by T test | No, only present if *de_analysis* is set
output_name.loom | Results in loom format | No, only present if *output_loom* is set
output_name_var.loom | Results containing only variable genes in loom format | No, only present if *output_loom* is set

