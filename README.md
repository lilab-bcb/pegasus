# scRNA-Seq

Table of Contents
-----------------

* [First Time Running - Authenticate with Google](#first-time-running---authenticate-with-google)
* [Run Cell Ranger mkfastq/count](#run-cell-ranger-mkfastqcount)
  * [Cell Ranger mkfastq/count inputs](#cell-ranger-mkfastqcount-inputs)
  * [Cell Ranger mkfastq/count outputs](#cell-ranger-mkfastqcount-outputs)
  * [Only run *cellranger count*](#only-run-cellranger-count)
* [Run Single Cell RNA-Seq analysis tools (scrtools)](#run-single-cell-rna-seq-analysis-tools-scrtools)
  * [scrtools steps](#scrtools-steps)
  * [global inputs](#global-inputs)
  * [aggregate_matrix](#aggregate_matrix)
  * [cluster](#cluster)
  * [de_analysis](#de_analysis)
  * [plot](#plot)
* [Run subcluster analysis](#run-subcluster-analysis)
  * [scrtools_subcluster steps](#scrtools_subcluster-steps)
  * [scrtools_subcluster's inputs](#scrtools_subclusters-inputs)
  * [scrtools_subcluster's outputs](#scrtools_subclusters-outputs)

## First Time Running - Authenticate with Google

1. Ensure the Google Cloud SDK is installed on your computer. 
    
    Note: Broad users do not have to install this-they can type *use Google-Cloud-SDK* to make the Google Cloud tools available. 

1. Type *gcloud auth login* to login to Google Cloud.

1. Copy and paste the link in your unix terminal into your web browser.

1. Enter authorization code in unix terminal


## Run Cell Ranger mkfastq/count
1. [Create](https://software.broadinstitute.org/firecloud/documentation/article?id=10746) a [FireCloud](https://portal.firecloud.org/) workspace or use an existing workspace.

1. Copy your sequencing output to the workspace bucket using gsutil in your unix terminal. 

    It is highly recommended that you delete the **BCL files** after the pipeline is finished by turning on the **delete_input_directory** option.

    Example of copying a directory to a Google Cloud bucket: *gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4*
    
    -m means copy in parallel, -r means copy the directory recursively.
    
    Note: Broad users need to be on an UGER node (not a login node) in order to use the -m flag
    
    You can also read [FireCloud instructions](https://software.broadinstitute.org/firecloud/documentation/article?id=10574) on uploading data.

1. Create a scRNA-Seq formatted sample sheet. Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.  
The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices.Note that *Sample*, *Lane*, and *Index* columns are defined exactly the same as in 10x's simple CSV layout file.

scRNA-Seq formatted sample sheet description (required column headers are shown in bold):

Column | Description
--- | --- | 
**Sample** | Contains sample names. Each 10x channel should have a unique sample name.
**Reference** | Provides the reference genome used by *cellranger count* for each 10x channel. The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as **GRCh38** for human, **mm10** for mouse, and **GRCh38\_and\_mm10** for human and mouse)
**Flowcell** | Indicates the Google bucket URL of uploaded BCL folders.
**Lane** | Tells which lanes the sample was pooled into.
**Index** | Contains 10x sample index set names (e.g. SI-GA-A12).
Chemistry | Optionally describe the 10x chemistry used for the sample. If this column is omitted, *cellranger count* will try to determine the chemistry automatically.

The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.
   
   Example:
   
    ```
    Sample,Reference,Flowcell,Lane,Index,Chemistry
    sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,threeprime
    sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime
    sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,threeprime
    sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime
    ```

1. Upload your sample sheet to the workspace bucket.  

    Example: *gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

1. Import cellranger_mkfastq_count method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type cellranger_mkfastq_count.

1. Uncheck "Configure inputs/outputs using the Workspace Data Model"

### Cell Ranger mkfastq/count inputs:

*Cell Ranger mkfastq/count* takes Illumina outputs as input and runs *cellranger mkfastq* and *cellranger count*. Please see the description of inputs below. Note that required inputs are shown in bold.

Name | Description | Example | Default
--- | --- | --- | ---
**input\_csv\_file** | Sample Sheet (contains Sample, Flowcell, Lane, Index) | "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv" | 
**cellranger\_output\_directory** | Cellranger output directory | "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output" |
run_mkfastq | If you want to run *cellranger mkfastq* | true | true
run_count | If you want to run *cellranger count* | true | true
delete\_input\_directory | If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges. | true | false
do_force_cells | force cells | true | true
force_cells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expect_cells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
cellranger_version | Cellranger version | "2.1.1" | "2.1.1"
num_cpu | Number of cpus to request for one node | 64 | 64
memory | Memory in GB | 128 | 128
mkfastq_disk_space | Optional disk space in gigabytes for mkfastq. | 1500 | 1500 
count_disk_space | Disk space in gigabytes needed for cell ranger count | 500 | 500
preemptible | Number of preemptible tries | 2 | 2

### Cell Ranger mkfastq/count outputs:

See the table below for important *Cell Ranger mkfastq/count* outputs.

Name | Type | Description
--- | --- | ---
output_fastqs_directory | Array[String] | A list of google bucket urls containing FASTQ files, one url per flowcell.
output_count_directory | Array[String] | A list of google bucket urls containing count matrices, one url per sample.
metrics_summaries | File | A excel spreadsheet containing QCs for each sample.
output_web_summary | Array[File] | A list of htmls visualizing QCs for each sample (cellranger count output).

### Only run *cellranger count*

Sometimes, people might want to perform demux locally and only run *cellranger count* on the cloud. This section describes how to perform count only use *cellranger_mkfastq_count*.

1. Copy your FASTQ files to the workspace using gsutil in your unix terminal. 

    You should upload folders of FASTQS. Each foloder contains one subfolder per sample. Each subfolder contains all FASTQ files from the corresponding sample.

    Example: *gsutil -m cp -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq*
    
    -m means copy in parallel, -r means copy the directory recursively.
    
    Note: Broad users need to be on an UGER node (not a login node) in order to use the -m flag
    
    You can also read [FireCloud instructions](https://software.broadinstitute.org/firecloud/documentation/article?id=10574) on uploading data.

1. Replace the *Flowcell* column in the sample sheet with the locations of FASTQ folders and upload it into your workspace.
   
   Example: 
    ```
    Lane,Sample,Index
    Sample,Reference,Flowcell,Lane,Index,Chemistry
    sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq,3-4,SI-GA-B8,threeprime
    sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq,7-8,SI-GA-D8,fiveprime
    sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/K10WBC9Z2_fastq,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/K10WBC9Z2_fastq,3-4,SI-GA-B8,threeprime
    sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/K10WBC9Z2_fastq,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/K10WBC9Z2_fastq,7-8,SI-GA-D8,fiveprime
    ```
1. Set optional input **run_mkfastq** to **false**.

## Run Single Cell RNA-Seq analysis tools (scrtools)

1. Create a sample sheet, **count_matrix.csv**, which describes the metadata for each 10x channel. The sample sheet should at least contain 3 columns --- *Sample*, *Reference*, and *Location*. *Sample* refers to sample names, *Reference* refers to the genome name, and *Location* refers to the location of the channel-specific count matrix in 10x format (e.g. gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5). You are free to add any other columns and these columns will be used in selecting channels for futher analysis. In the example below, we have *Source*, which refers to the tissue of origin, *Platform*, which refers to the sequencing platform, and *Donor*, which refers to the donor ID.

    Example:
    ```
    Sample,Source,Platform,Donor,Reference,Location
    sample_1,bone_marrow,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5
    sample_2,bone_marrow,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/filtered_gene_bc_matrices_h5.h5
    sample_3,pbmc,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_3/filtered_gene_bc_matrices_h5.h5
    sample_4,pbmc,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_4/filtered_gene_bc_matrices_h5.h5
    ```    
1. Upload your sample sheet to the workspace.  

    Example: *gsutil cp /foo/bar/projects/my_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

1. Import **scrtools** method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type **scrtools**.

1. Uncheck "Configure inputs/outputs using the Workspace Data Model"

### scrtools steps:

**scrtools** processes single cell data in the following steps:

1. **aggregate_matrix**. This step aggregates channel-specific count matrices into one big count matrix. Users could specify which channels they want to analyze and which sample attributes they want to import to the count matrix in this step.

1. **cluster**. This step is the main analysis step. In this step, **scrtools** performs low quality cell filtration, variable gene selection, batch correction, dimension reduction, diffusion map calculation, graph-based clustering and 2D visualization calculation (e.g. tSNE/FLE).

1. **de_analysis**. This step is optional. In this step, **scrtools** could calculate potential markers for each cluster by performing a variety of differential expression (DE) analysis. The available DE tests include Welch's t test, Fisher's exact test, and Mann-Whitney U test. **scrtools** could also calculate the area under ROC curve values for putative markers. If the samples are human or mouse immune cells, **scrtools** could also optionally annotate putative cell types for each cluster based on known markers.

1. **plot**. This step is optional. In this step, **scrtools** could generate 3 types of figures based on the **cluster** step results. First, **composition** plots are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects. Second, **tsne** plot shows the same tSNE colored by different attributes (e.g. cluster labels, conditions) side-by-side. Lastly, **diffmap** plots are 3D interactive plots showing the diffusion maps. The 3 coordinates are the first 3 PCs of all diffusion components.

In the following, we will first introduce global inputs and then introduce the WDL inputs and outputs for each step separately. But please note that you need to set inputs from all steps simultaneously in the FireCloud WDL. 

Note that we will make the required inputs/outputs bold and all other inputs/outputs are optional.


### global inputs

Name | Description | Example | Default
--- | --- | --- | ---
**input_count_matrix_csv** | Input CSV file describing metadata of each 10x channel | "my_count_matrix.csv" | 
**output_name** | This is the prefix for all output files. It should contain the google bucket url, subdirectory name and output name prefix | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results" | 
genome | Reference genome name | "mm10" | "GRCh38"
num_cpu | Number of cpus per scrtools job | 32 | 64
memory | Memory size in GB | 200 | 200
disk_space | Total disk space | 100 | 100
preemptible | Number of preemptible tries | 2 | 2

### aggregate_matrix

#### aggregate_matrix inputs

Name | Description | Example | Default
--- | --- | --- | ---
restrictions | Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. Multiple restrictions are separated by ';' | "Source:bone_marrow;Platform:NextSeq" | 
attributes | Specify a comma-separated list of outputted attributes. These attributes should be column names in the count_matrix.csv file | "Source,Platform,Donor" | 

#### aggregate_matrix output

Name | Type | Description
--- | --- | ---
**output_10x_h5** | File | Aggregated count matrix in 10x format

### cluster

#### cluster inputs

Note that we will only list important inputs here. For other inputs, please refer to **scrtools** package documentation.

Name | Description | Example | Default
--- | --- | --- | ---
output_filtration_results | If output cell and gene filtration results to a spreadsheet | true | true
output_loom | If output loom-formatted file | false | false
correct_batch_effect | If correct batch effects | false | false
batch_group_by | Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either ‘attr’, or ‘attr1+attr2+…+attrn’, or ‘attr=value11,…,value1n_1;value21,…,value2n_2;…;valuem1,…,valuemn_m’. In the first form, ‘attr’ should be an existing sample attribute, and groups are defined by ‘attr’. In the second form, ‘attr1’,…,’attrn’ are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute ‘attr’ has a value among valuei1,…,valuein_i. A cell belongs to group 0 if it does not belong to any other groups | "Donor" | None
min_genes | Only keep cells with at least <number> of genes | 500 | 500
max_genes | Only keep cells with less than <number> of genes | 6000 | 6000
mito_prefix | Prefix for mitochondrial genes | "mt-" | "MT-"
percent_mito | Only keep cells with mitochondrial ratio less than <ratio> | 0.1 | 0.1 
gene_percent_cells | Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes | 0.0005 | 0.0005
counts_per_cell_after | Total counts per cell after normalization | 1e5 | 1e5
random_state | Random number generator seed | 0 | 0
nPC | Number of principal components | 50 | 50
nDC | Number of diffusion components | 50 | 50
diffmap_K | Number of neighbors used for constructing affinity matrix | 100 | 100
diffmap_alpha | Power parameter for diffusion-based pseudotime | 0.5 | 0.5 
run_louvain | Run louvain clustering algorithm | true | false
louvain_resolution | Resolution parameter for the louvain clustering algorithm | 1.3 | 1.3
run_approximated_louvain | Run approximated louvain clustering algorithm | true | false
approx_louvain_ninit | Number of Kmeans tries | 30 | 20 
approx_louvain_nclusters | Number of clusters for Kmeans initialization | 40 | 30
approx_louvain_resolution | Resolution parameter for louvain | 1.3 | 1.3
run_tsne | Run multi-core tSNE for visualization | true | false
tsne_perplexity | tSNE’s perplexity parameter | 30 | 30
run_fitsne | Run FItSNE for visualization | true | false
run_umap | Run umap for visualization | true | false
umap_on_diffmap | Run umap on diffusion components | "ture" | false
run_fle | Run force-directed layout embedding | true | false
fle_K | K neighbors for building graph for FLE | 50 | 50
fle_n_steps | Number of iterations for FLE | 10000 | 10000

#### cluster outputs

Name | Type | Description
--- | --- | ---
**output_h5ad** | File | h5ad-formatted HDF5 file containing all results (output_name.h5ad)
output_filt_xlsx | File | Spreadsheet containing filtration results (output_name.filt.xlsx)
output_loom_file | File | Outputted loom file (output_name.loom)

### de_analysis

#### de_analysis inputs

Name | Description | Example | Default
--- | --- | --- | ---
perform_de_analysis | If perform de analysis | true | true
cluster_labels | Specify the cluster labels used for differential expression analysis | "louvain_labels" | "louvain_labels" 
alpha | Control false discovery rate at <alpha> | 0.05 | 0.05
fisher | Calculate Fisher’s exact test | true | false
mwu | Calculate Mann-Whitney U test | true | false
roc | Calculate area under cuver in ROC curve | true | false
annotate_cluster | If also annotate cell types for clusters based on DE results | true | false
organism | Organism, could either be "human" or "mouse" | "mouse" | "human"
minimum_report_score | Minimum cell type score to report a potential cell type | 0.5 | 0.5

#### de_analysis outputs

Name | Type | Description
--- | --- | ---
output_de_h5ad | File | h5ad-formatted results with DE results updated (output_name.h5ad)
output_de_xlsx | File | Spreadsheet reporting DE results (output_name.de.xlsx)
output_anno_file | File | Annotation file (output_name.anno.txt)

### plot

#### plot inputs

Name | Description | Example | Default
--- | --- | --- | ---
plot_composition | Takes the format of "label:attr,label:attr,...,label:attr". If non-empty, generate composition plot for each "label:attr" pair. "label" refers to cluster labels and "attr" refers to sample conditions | "louvain_labels:Donor" | None
plot_tsne | Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored tSNEs side by side | "louvain_labels,Donor" | None
plot_diffmap | Takes the format of "attr,attr,...,attr". If non-empty, generate attr colored 3D interactive plot. The 3 coordinates are the first 3 PCs of all diffusion components | "louvain_labels,Donor" | None

#### plot outputs

Name | Type | Description
--- | --- | ---
output_pngs | Array[File] | Outputted png files
output_htmls | Array[File] | Outputted html files


## Run subcluster analysis

Once we have **scrtools** outputs, we could further analyze a subset of cells by running **scrtools_subcluster**. To run **scrtools_subcluster**, follow the following steps:

1. Import **scrtools_subcluster** method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type **scrtools_subcluster**.

1. Uncheck "Configure inputs/outputs using the Workspace Data Model"

### scrtools_subcluster steps:

*scrtools_subcluster* processes the subset of single cells in the following steps:

1. **subcluster**. In this step, **scrtools_subcluster** first select the subset of cells from **scrtools** outputs according to user-provided criteria. It then performs batch correction, dimension reduction, diffusion map calculation, graph-based clustering and 2D visualization calculation (e.g. tSNE/FLE).

1. **de_analysis**. This step is optional. In this step, **scrtools_subcluster** could calculate potential markers for each cluster by performing a variety of differential expression (DE) analysis. The available DE tests include Welch's t test, Fisher's exact test, and Mann-Whitney U test. **scrtools_subcluster** could also calculate the area under ROC curve values for putative markers. If the samples are human or mouse immune cells, **scrtools_subcluster** could also optionally annotate putative cell types for each cluster based on known markers.

1. **plot**. This step is optional. In this step, **scrtools_subcluster** could generate 3 types of figures based on the **subcluster** step results. First, **composition** plots are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects. Second, **tsne** plot shows the same tSNE colored by different attributes (e.g. cluster labels, conditions) side-by-side. Lastly, **diffmap** plots are 3D interactive plots showing the diffusion maps. The 3 coordinates are the first 3 PCs of all diffusion components.

### scrtools_subcluster's inputs

Since **scrtools_subcluster** shares many inputs/outputs with **scrtools**, we will only cover inputs/outputs that are specific to **scrtools_subcluster**.

Note that we will make the required inputs/outputs bold and all other inputs/outputs are optional.

Name | Description | Example | Default
--- | --- | --- | ---
**input_h5ad** | Input h5ad file containing *scrtools* results | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results.h5ad" | 
**output_name** | This is the prefix for all output files. It should contain the google bucket url, subdirectory name and output name prefix | "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results_sub" | 
**subset_selections** | Specify which cells will be included in the subcluster analysis. This field contains one or more <subset_selection> strings separated by ';'. Each <subset_selection> string takes the format of ‘attr:value,…,value’, which means select cells with attr in the values. If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings | "louvain_labels:3,6" | 
calculate_pseudotime | Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes | "sample_1-ACCCGGGTTT-1" | None
num_cpu | Number of cpus per scrtools job | 32 | 64
memory | Memory size in GB | 200 | 200
disk_space | Total disk space | 100 | 100
preemptible | Number of preemptible tries | 2 | 2

### scrtools_subcluster's outputs

Name | Type | Description
--- | --- | --- 
**output_h5ad** | File | h5ad-formatted HDF5 file containing all results (output_name.h5ad)
output_loom_file | File | Outputted loom file (output_name.loom)
output_de_h5ad | File | h5ad-formatted results with DE results updated (output_name.h5ad)
output_de_xlsx | File | Spreadsheet reporting DE results (output_name.de.xlsx)
output_pngs | Array[File] | Outputted png files
output_htmls | Array[File] | Outputted html files

