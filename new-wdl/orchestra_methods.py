import os
from subprocess import call
import xml.etree.ElementTree as ET
import pandas as pd
from subprocess import run
from subprocess import check_output
import sys
import argparse



def cellranger_count(sampleID, commaFastqs='',fastqs='', expectCells='', forceCells='', secondary='true', chemistry='threeprime'):
    dirs = set()
    sample = sampleID
    i = 0
    if commaFastqs is not '':
        fastqs = commaFastqs.split(",")
    for f in fastqs:
        # get the fastqs
        call(["mkdir", str(i)])
        call(["gsutil", "-q", "-m", "cp", "-r", f, str(i)])
        os.path.join(str(i), sample,"")
        dirs.add(os.path.join(str(i), sample,""))
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

def cellranger_mkfastq(bcl,masterCsv,output_directory):
    # open up log.txt
    df = pd.read_csv('log.txt', header=0)
    path = bcl
    run = list(filter(None, path.split("/")))[-1]

    # get flowcell
    tree = ET.parse(os.path.join(run,"RunInfo.xml"))
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
    call_args.append('--run=' + os.path.join(run,""))
    call_args.append('--csv=' + 'sample.csv')
    call_args.append('--output-dir=fastqs')
    call(call_args)

    # move the fastqs to the output directory
    call(["gsutil", "-q", "-m", "mv", os.path.join('fastqs',flowcell), output_directory])
    call(["gsutil", "-q", "-m", "mv", os.path.join(flowcell,"outs","qc_summary.json"), os.path.join(output_directory,flowcell+"_qc_summary.json")])
    
    file = open("path.txt","w") 
    file.write(os.path.join(output_directory,flowcell,""))
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
      key = os.path.join(sample,"")
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

def __main__(argv):
    command = argv[1].replace("-c=", "")
    parser = argparse.ArgumentParser()
    if command == "count":
        parser.add_argument('--sampleId', '-id', help="Id of sample being run", type= str)
        parser.add_argument('--commaFastqs', '-cf', help="Comma seperated String with list of fastq directories", type= str, default='')
        parser.add_argument('--fastqs', '-fs', help="List of fastq directories", type= list, default=[''])
        parser.add_argument('--expectCells', '-E', help="Number of cells to expect", type= str, default='')
        parser.add_argument('--forceCells', '-F', help="Force number of cells", type= str, default='')
        parser.add_argument('--chemistry', '-C', help="Chemistry of fastqs", type= str, default = "threeprime")
        parser.add_argument('--secondary', '-S', help="Run cellranger secondary analysis", type= str, default = "true")
        parser.add_argument('--command', '-c', help="Command to run", type= str, default = "")

        args = parser.parse_args()
        cellranger_count(sampleId= args.sampleId, commaFastqs= args.commaFastqs, fastqs= args.fastqs, expectCells= args.expectCells, forceCells= args.forceCells, chemistry = args.chemistry)
    elif command == "mkfastq":
        parser = argparse.ArgumentParser()
        parser.add_argument('--bcl', '-b', help="Location of bcl", type= str)
        parser.add_argument('--masterCsv', '-M', help="Master Csv file containing maps of information", type= str, default='')
        parser.add_argument('--output_directory', '-O', help="List of fastq directories", type= str, default='')
        parser.add_argument('--command', '-c', help="Command to run", type= str, default = "")
        
        args = parser.parse_args()
        cellranger_mkfastq(bcl = args.bcl, masterCsv = args.masterCsv, output_directory = args.output_directory)
    elif command == "parse":
        parser = argparse.ArgumentParser()
        parser.add_argument('--masterCsv', '-M', help="Master Csv file containing maps of information", type= str, default='')
        parser.add_argument('--command', '-c', help="Command to run", type= str, default = "")
        
        args = parser.parse_args()
        orchestra_parse_csv(masterCsv = args.masterCsv)
    elif command == "analysis":
        parser = argparse.ArgumentParser()
        parser.add_argument('--masterCsv', '-M', help="Master Csv file containing maps of information", type= str, default='')
        parser.add_argument('--h5s', '-hs', help="H5 output files", type= list, default=[''])
        parser.add_argument('--command', '-c', help="Command to run", type= str, default = "")
        
        args = parser.parse_args()
        orchestra_analysis_csv(masterCsv = args.masterCsv, h5s = args.h5s)
    elif command == "filter":
        parser = argparse.ArgumentParser()
        parser.add_argument('--masterCsv', '-M', help="Master Csv file containing maps of information", type= str, default='')
        parser.add_argument('--paths', '-p', help="Paths to fastq directories", type= list, default=[''])
        parser.add_argument('--sampleIds', '-S', help="List of Sample Names", type= list, default=[''])
        parser.add_argument('--command', '-c', help="Command to run", type= str, default = "")
        
        args = parser.parse_args()
        orchestra_filter(masterCsv = args.masterCsv, paths = args.paths, sampleIds = args.sampleIds)
    else:
        print("Error", command, "Is not a registered command")
if __name__ == "__main__":
    __main__(sys.argv)