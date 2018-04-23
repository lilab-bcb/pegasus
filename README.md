# scRNA-Seq

## First Time Running - Authenticate with Google

1. Ensure the Google Cloud SDK is installed. 
    
    Note: Broad users can type *use Google-Cloud-SDK* to make the Google Cloud tools available. 

1. Type *gcloud auth login*

1. Open link in browser

1. Enter authorization code in unix terminal


## Run Cell Ranger mkfastq/count
1. Create a FireCloud workspace or use an existing workspace.

1. Upload your sequencing output to the workspace.

    Example: *gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/VK18WBC6Z4*

1. Upload your CSV file to the workspace

    Example: *gsutil cp /foo/bar/projects/VK18WBC6Z4/my_bcl.csv gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/*

1. Upload "fake" data module required by FireCloud

    Download https://github.com/broadinstitute/scRNA-Seq/raw/master/wdl/participant_model.txt

    In FireCloud select the "Data" tab. Click "Import Metadata". Select "Import From File". Select the file and click "Upload"
1. Import cellranger_mkfastq_count method.
    
    In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Copy from another workspace". Type cellranger_mkfastq_count.

### Cell Ranger mkfastq/count inputs:

Name | Description | Example | Default
--- | --- | --- | ---
input_csv_file | 3 column CSV (Lane, Sample, Index) | "gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/my_file.csv"
input_directory | Sequencer output directory containing Config/, Data/, Images/, InterOp/, etc. | "gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/my_dir" | 
fastq_output_directory | Fastq output directory | "gs://fc-e7760257-20a8-46ae-8358-93a83f02a0ba/my_dir" |
mkfastq_disk_space | Optional disk space for mkfastq. | 500 | 500 
cell_ranger_version | Cell ranger version | "2.1.1" | "2.1.1"
transcriptome | mm10, GRCh38, or a gs URL to a transcriptome directory tar.gz | "mm10" | 
secondary | Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.) | false | false
do_force_cells | force cells | true | true
force_cells | Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells | 3000 | 6000
expect_cells | Expected number of recovered cells. Mutually exclusive with force_cells | 1000 | 3000
count_disk_space | Disk space needed for cell ranger count | 500 | 250

