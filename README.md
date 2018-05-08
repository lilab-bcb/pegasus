# scRNA-Seq

Table of Contents
-----------------

* [First Time Running - Authenticate with Google](#first_time)
* [Run Cell Ranger mkfastq/count](#run_cellranger)
  * [Cell Ranger mkfastq/count inputs](#cellranger_mkfastq_count)
* [Run Single Cell RNA-Seq analysis tools](#run_scrtools)
  * [Prepare count_matrix.csv](#count_matrix_csv)
  * [scrtools_aggregate_matrix](#aggr_mat)
  * [scrtools_merge_matrix](#merge_mat)
  * [scrtools_cluster](#cluster)
  * [scrtools_annotate](#annotate)
  * [scrtools_subcluster](#subcluster)

## <a name="first_time"></a> First Time Running - Authenticate with Google

1. Ensure the Google Cloud SDK is installed. 
    
    Note: Broad users can type *use Google-Cloud-SDK* to make the Google Cloud tools available. 

1. Type *gcloud auth login*

1. Open link in browser

1. Enter authorization code in unix terminal


## <a name="run_cellranger"></a> Run Cell Ranger mkfastq/count
1. Create a FireCloud workspace or use an existing workspace.

1. Upload your sequencing output to the workspace.

    Example: *gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4*

1. Upload your BCL CSV file to the workspace

    Example: *gsutil cp /foo/bar/projects/my_bcl.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

1. Upload "dummy" data model required by FireCloud

    Download https://github.com/broadinstitute/scRNA-Seq/raw/master/wdl/participant_model.txt

    In FireCloud select the "Data" tab. Click "Import Metadata". Select "Import From File". Select the file and click "Upload"
1. Import cellranger_mkfastq_count method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Copy from another workspace". Type cellranger_mkfastq_count.

### <a name="cellranger_mkfastq_count"></a> Cell Ranger mkfastq/count inputs:

*Cell Ranger mkfastq/count* takes Illumina outputs as its input and run *cellranger mkfastq* and *cellranger count* for you. Please see the description of inputs below. Note that required inputs are shown in bold.

Name | Description | Example | Default
--- | --- | --- | ---
**input_csv_file** | 3 column CSV (Lane, Sample, Index) | "gs://fc-e0000000-0000-0000-0000-000000000000/my_file.csv" | 
**input_directory** | Sequencer output directory containing Config/, Data/, Images/, InterOp/, etc. | "gs://fc-e0000000-0000-0000-0000-000000000000/my_dir" | 
**cellranger_output_directory** | Cellranger output directory | "gs://fc-e0000000-0000-0000-0000-000000000000/my_dir" |
**transcriptome** | gs URL to a transcriptome tar.gz ("gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz" for mm10, "gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz" for GRCh38 | "gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz" | 
mkfastq_disk_space | Optional disk space for mkfastq. | 500 | 500 
cellranger_version | Cellranger version | "2.1.1" | "2.1.1"
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
do_force_cells | force cells | true | true
force_cells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expect_cells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
count_disk_space | Disk space needed for cell ranger count | 500 | 250


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

In the above CSV file, *Sample*, *Reference*, and *Location* are required. *Sample* refers to the sample names listed in the BCL CSV file, *Reference* refers to the genome name, and *Location* refers to the cellranger output folder. The cellranger count outputs are expected to locate at *Location/Sample* (e.g. gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S1). You are free to add any other columns and these columns will be used in selecting channels for futher analysis. In this example, we have *Source*, which refers to the tissue of origin, *Platform*, which refers to the sequencing platform, and *Donor*, which refers to donor ID.

You should upload **count_matrix.csv** to your workspace

    Example: *gsutil cp /foo/bar/projects/my_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

Then you can aggregate 10x count matrices into a single count matrix using **scrtools_aggregate_matrix**

### <a name="aggr_mat"></a> scrtools_aggregate_matrix

*scrtools_aggregate_matrix* is used to aggregate individual 10x channels into a big count matrix for downstream analysis. Please see the inputs below.

Name | Description | Example | Default
--- | --- | --- | ---
**input_count_matrix_csv** | Input CSV file describing metadata of each 10x channel | "my_count_matrix.csv" | 
**output_folder** | This is the folder for all analysis results | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**output_name** | Output file name of this task, the count matrix *output_name_10x.h5* and metadata file *output_name.attr.csv* will be generated | "my_aggr_mat_bm" | 
genome | The genome cellranger used to generate count matrices | "GRCh38" | "GRCh38"
restrictions | Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. Multiple restrictions are separated by ';'. If not restrictions are provided, all channels in the *count_matrix.csv* will be selected | "Source:bone_marrow" | 
attributes | Specify a comma-separated list of outputted attributes (i.e. column names in *count_matrix.csv*). These attributes will be imported to *output_name.attr.csv* | "Source,Donor" | 
groupby | When we know there are different groups in the study, such as bone_marrow and pbmc, we could perform batch correction separately for each group. This optional field is used to create a new attribute, GroupBy, that helps to identify groups. It takes the format of 'value' or 'attribute=value' or 'attr1+attr2+....' | "Source+Donor" | 
diskSpace | Disk space needed for this task | 100 | 100

### <a name="merge_mat"></a> scrtools_merge_matrix

*scrtools_merge_matrix* is used to combine two or more aggregated matrices produced by *scrtools_aggregate_matrix* into one big count matrix. This step is optional. 

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools_aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**input_names** | A comma-separated list of input names for count matrices that you want to merge | "my_aggr_mat_bm,my_aggr_mat_pbmc" | 
**output_name** | Output file name of this task, the count matrix *output_name_10x.h5* and metadata file *output_name.attr.csv* will be generated | "my_aggr_mat" | 
genome | The genome cellranger used to generate count matrices | "GRCh38" | "GRCh38"
symbols | A comma-separated list of symbols representing each input matrix, this input is used with 'attributes' | "bm,pbmc" | 
attributes | A comma-separated list of attributes. When merging matrices, the matrix symbol defined in 'symbols' will be added in front of these attributes | "Donor" | 
diskSpace | Disk space needed for this task | 100 | 100

### <a name="cluster"></a> scrtools_cluster

*scrtools_cluster* performs PCA, tSNE visualization, diffusion map, Louvain clustering, and differential expression analysis (fisher's exact test and t test). This is the main step people want to know. Please see inputs below.

Note: currently, this step uses Scanpy as its backend engine. In the future, we will use more optimized codes to replace Scanpy.

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools_aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**input_name** | Input name of the aggreagated/merged matric | "my_aggr_mat" | 
**output_name** | Output file name of this task | "results" | 
num_cpu | Number of CPUs to use | 64 | 64
genome | The genome cellranger used to generate count matrices | "GRCh38" | "GRCh38"
output_filtration_results | A boolean value indicating if you want to output QC metrics into a spreadsheet | "true" | "false"
correct_batch_effect | A boolean value indicating if you want to correct for batch effects | "true" | "false"
batch_group_by | If correct for batch effects, if you want to correct it separately for each group. Group is defined by this input. It should be an attribute (e.g. GroupBy) | "Source" | 
plot_by_side | By default, this step will generate one tSNE plot colored by louvain cluster labels. If you provide an attribute name here, it will put the same tSNE colored by the attribute on the right-hand side | "Source" | 
legend_on_data | A boolean value indicating if we should put legends on data | "true" | "false"
plot_composition | A boolean value indicating if you want to generate composition plots for the attribute in *plot_by_side*, which for each cluster shows the percentage of cells in each attribute value | "true" | "false"
figure_size | Sizes of the composition plots in inches | "10,8" | "6,4"
plot_diffusion_map | A boolean value indicating if interactive 3D diffusion maps should be generated | "true" | "false" 
de_analysis | A boolean value indicating if you want to perform differential expression analysis | "true" | "true"
fold_change | Minimum fold change in either percentage (fisher test) or log expression (t test) to report a DE gene in spreadsheets | 2 | 1.5
labels | Cluster labels that help de_analysis identify clusters | "louvain_labels" | "louvain_labels"
output_loom | A boolean value indicating if you want to output loom-formatted files | "false" | "false"
import_attributes | Import attributes contained in the comma-separated list into the analysis object | "Donor" | 
min_genes | The minimum number of expressed genes to be consider a valid cell | 500 | 500
max_genes | *max_genes* - 1 is the maximum number of expressed genes to be consider a valid cell | 6000 | 6000
mito_prefix | Prefix of mitochondrial genes | "MT-" | "MT-"
percent_mito | Only keep cells with mitochondrial ratio less than *percent_mito* | 0.1 | 0.1
gene_percent_cells | Genes expressed in with less than *gene_percent_cells* * *number_of_cells* cells will be excluded from variable gene selection step. The default value requires a gene expressed in at least 3 cells when there are 6,000 cells | 0.0005 | 0.0005
counts_per_cell_after | Normalize each cell so that the sum of its normalized counts is equal to *counts_per_cell_after* | 1e5 | 1e5
louvain_resolution | Resolution parameter of the louvain clustering algorithm | 1 | 1.3
diskSpace | Disk space needed for this task | 250 | 250

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

### <a name="annotate"></a> scrtools_annotate

This step outputs putative cell type annotations for each cluster based on known markers. It required differential expression analyses performed in *scrtools_cluster* step. This step is optional and currently it only works for human immune cells. Please see the inputs below.

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools_aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**file_name** | Input file name, *file_name_de.h5ad* should exist | "results" | 
minimum_report_score | This step calculates a score between [0, 1] for each pair of cluster and putative cell type. It only report putative cell types with a minimum score of *minimum_report_score* | 0.1 | 0.5
no_use_non_de | Do not consider non-differentially-expressed genes as down-regulated | "true" | "false"
diskSpace | Disk space needed for this task | 100 | 100

### <a name="subcluster"></a> scrtools_subcluster

This step performs subcluster analysis based on *scrtools_cluster* outputs. Please see the inputs below.

Name | Description | Example | Default
--- | --- | --- | ---
**data_folder** | This is the folder for all analysis results. It should be the same one used in *scrtools_aggregate_matrix* | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir" | 
**input_name** | Input name of *scrtools_cluster* result | "results" | 
**output_name** | Output file name of this task | "results_sub" | 
**cluster_ids** | A comma-separated list of cluster IDs (numbered from 1) for subcluster | "1,3,9" | 
num_cpu | Number of CPUs to use | 64 | 64
correct_batch_effect | A boolean value indicating if you want to correct for batch effects | "true" | "false"
plot_by_side | By default, this step will generate one tSNE plot colored by louvain cluster labels. If you provide an attribute name here, it will put the same tSNE colored by the attribute on the right-hand side | "Source" | 
legend_on_data | A boolean value indicating if we should put legends on data | "true" | "false"
plot_composition | A boolean value indicating if you want to generate composition plots for the attribute in *plot_by_side*, which for each cluster shows the percentage of cells in each attribute value | "true" | "false"
figure_size | Sizes of the composition plots in inches | "10,8" | "6,4"
plot_diffusion_map | A boolean value indicating if interactive 3D diffusion maps should be generated | "true" | "false" 
de_analysis | A boolean value indicating if you want to perform differential expression analysis | "true" | "true"
fold_change | Minimum fold change in either percentage (fisher test) or log expression (t test) to report a DE gene in spreadsheets | 2 | 1.5
labels | Cluster labels that help de_analysis identify clusters | "louvain_labels" | "louvain_labels"
louvain_resolution | Resolution parameter of the louvain clustering algorithm | 1 | 1.3
output_loom | A boolean value indicating if you want to output loom-formatted files | "false" | "false"
diskSpace | Disk space needed for this task | 250 | 250

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

