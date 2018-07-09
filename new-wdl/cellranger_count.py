import os
from subprocess import call
import xml.etree.ElementTree as ET
import pandas as pd
from subprocess import run
from subprocess import check_output

def cellranger_count(sampleID, commaFastqs='',fastqs='', expectCells, forceCells, secondary, chemistry):
    dirs = set()
    sample = sampleID
    i = 0
    if commaFastqs is not '':
        fastqs = commaFastqs.split(",")
    else:
        fastqs = fastqs
    for f in fastqs:
        # get the fastqs
        call(["mkdir", str(i)])
        call(["gsutil", "-q", "-m", "cp", "-r", f, str(i)])
        dirs.add(str(i)+"/" +sample+"/")
        i+=1
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

def cellranger_mkfastq(bcl,masterCsv,output_directory,):
    # open up log.txt
    df = pd.read_csv('log.txt', header=0)
    path = bcl
    run = list(filter(None, path.split("/")))[-1]

    # get flowcell
    tree = ET.parse(run+"/RunInfo.xml")
    root = tree.getroot()
    flowcell = root.find('Run').find('Flowcell').text

    # Create the csv
    df = pd.read_csv(masterCsv,header=0)
    df = df.loc[df['Flowcell'] == path]
    df = df[["Lane","Sample", "Index"]]
    df.to_csv('sample.csv',index=False)

    # run mkfastq command
    call_args = list()
    call_args.append('cellranger')
    call_args.append('mkfastq')
    call_args.append('--run=' + run+"/")
    call_args.append('--csv=' + 'sample.csv')
    call_args.append('--output-dir=fastqs')
    call(call_args)

    # move the fastqs to the output directory
    call(["gsutil", "-q", "-m", "mv", 'fastqs/'+flowcell, output_directory])
    call(["gsutil", "-q", "-m", "mv", flowcell+"/outs/qc_summary.json", output_directory+flowcell+"_qc_summary.json"])
    
    file = open("path.txt","w") 
    file.write(output_directory+flowcell+"/")
    file.close()

def orchestra_parse_csv(masterCsv):
    df = pd.read_csv(masterCsv,header=0)
    bcl_paths = set(df['Flowcell'].tolist())
    bcl_file = open('bcls.txt', 'w+')
    for item in bcl_paths:
      bcl_file.write("%s\n" % item)
    bcl_file.close()
    sampleIds = set(df['Sample'].tolist())
    samples_file = open('samples.txt', 'w+')
    for item in sampleIds:
      samples_file.write("%s\n" % item)
    samples_file.close()


def orchestra_analysis_csv(masterCsv, h5s):
    df = pd.read_csv(masterCsv,header=0)
    df = df.drop(columns = ["Flowcell", "Lane", "Index"])
    h5s = sorted(h5s)
    sampleIds = set(df['Sample'].tolist())
    #TODO Parse h5s to make sure it's matching the right one-- I think this is a good way to do it
    df = df.sort_values(by=["Sample"]).drop_duplicates(subset=["Sample"])
    # This is so ugly
    h5s = [h5.replace("/cromwell_root/", "gs://") for h5 in h5s]
    df["Location"] = h5s
    df.to_csv("analysis.csv", index=False)


def orchestra_filter(paths, masterCsv, sampleIds):
    fastqs = []
    for f in paths:
      fastqs = fastqs + check_output(["gsutil", "ls", f]).decode("utf-8").split("\n")

    chemistry_list = open('chemistry.tsv', 'w+')
    genome_list = open('genome.tsv', 'w+')
    df = pd.read_csv(masterCsv)
    paths_list = open('paths.tsv', 'w+')
    # Now we have a list of every sample
    for sample in sampleIds:
      # Sample Paths Map
      key = sample + "/"
      filter_paths = ",".join([path for path in fastqs if path.endswith(key)])
      paths_list.write("%s\t%s\n" % (sample, filter_paths))
      # Chemistry and Genome Map
      rows = df.loc[df["Sample"] == sample]
      chemistry = rows["Chemistry"].tolist()[0]
      genome = rows["Reference"].tolist()[0]
      # Write to files
      chemistry_list.write("%s\t%s\n" % (sample, chemistry))
      genome_list.write("%s\t%s\n" % (sample, genome))
    genome_list.close()
    chemistry_list.close()
    paths_list.close()