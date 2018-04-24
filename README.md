# scRNA-Seq

Table of Contents
-----------------

* [First Time Running - Authenticate with Google](#first_time)
* [Run Cell Ranger mkfastq/count](#run_cellranger)
* [Run Single Cell RNA-Seq analysis tools](#run_scrtools)

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

### Cell Ranger mkfastq/count inputs:

Name | Description | Example | Default
--- | --- | --- | ---
input_csv_file | 3 column CSV (Lane, Sample, Index) | "gs://fc-e0000000-0000-0000-0000-000000000000/my_file.csv"
input_directory | Sequencer output directory containing Config/, Data/, Images/, InterOp/, etc. | "gs://fc-e0000000-0000-0000-0000-000000000000/my_dir" | 
cellranger_output_directory | Cellranger output directory | "gs://fc-e0000000-0000-0000-0000-000000000000/my_dir" |
transcriptome | gs URL to a transcriptome tar.gz ("gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz" for mm10, "gs://regev-lab/resources/cellranger/refdata-cellranger-GRCh38-1.2.0.tar.gz" for GRCh38 | "gs://regev-lab/resources/cellranger/refdata-cellranger-mm10-1.2.0.tar.gz" | 
mkfastq_disk_space | Optional disk space for mkfastq. | 500 | 500 
cellranger_version | Cellranger version | "2.1.1" | "2.1.1"
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
do_force_cells | force cells | true | true
force_cells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expect_cells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
count_disk_space | Disk space needed for cell ranger count | 500 | 250

## <a name="run_scrtools"></a> Run Single Cell RNA-Seq analysis tools (scrtools)

Before we run the scrtools, we need to first prepare a CSV file, *count_matrix.csv*, which describes the metadata for each 10x channel.

### count_matrix.csv

```
Sample,Source,Platform,Donor,Reference,Location
S1,bone_marrow,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
S2,bone_marrow,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
S3,pbmc,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
S4,pbmc,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir
```

In the above CSV file, fileds *Sample*, *Reference*, and *Location* are required. *Sample* refers to the sample names listed in the BCL CSV file, *Reference* refers to the genome name, and *Location* refers to the cellranger output folder. The cellranger count outputs are expected to locate at *Location/Sample* (e.g. gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/S1). You are free to add any other columns and these columns will be used in selecting channels for futher analysis. In this example, we have *Source*, which refers to the tissue of origin, *Platform*, which refers to the sequencing platform, and *Donor*, which refers to donor ID.

You should upload **count_matrix.csv** to your workspace

    Example: *gsutil cp /foo/bar/projects/my_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/*

Then you can aggregate 10x count matrices into a single count matrix using **scrtools_merge_matrix**

### scrtools_merge_matrix

