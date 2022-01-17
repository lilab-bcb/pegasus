"""
Pegasus is a tool for analyzing transcriptomes of millions of single cells. It is a command line tool, a python package and a base for Cloud-based analysis workflows.

Usage:
  pegasus <command> [<args>...]
  pegasus -h | --help
  pegasus -v | --version

Sub-commands:
  Preprocessing:
    aggregate_matrix        Aggregate cellranger-outputted channel-specific count matrices into a single count matrix. It also enables users to import metadata into the count matrix.
  Demultiplexing:
    demuxEM                 Demultiplex cells/nuclei based on DNA barcodes for cell-hashing and nuclei-hashing data.
  Analyzing:
    cluster                 Perform first-pass analysis using the count matrix generated from 'aggregate_matrix'. This subcommand could perform low quality cell filtration, batch correction, variable gene selection, dimension reduction, diffusion map calculation, graph-based clustering, visualization. The final results will be written into zarr-formatted file, or h5ad file, which Seurat could load.
    de_analysis             Detect markers for each cluster by performing differential expression analysis per cluster (within cluster vs. outside cluster). DE tests include Welch's t-test, Fisher's exact test, Mann-Whitney U test. It can also calculate AUROC values for each gene.
    find_markers            Find markers for each cluster by training classifiers using LightGBM.
    annotate_cluster        This subcommand is used to automatically annotate cell types for each cluster based on existing markers. Currently, it works for human/mouse immune/brain cells, etc.
  Plotting:
    plot                    Make static plots, which includes plotting tSNE, UMAP, and FLE embeddings by cluster labels and different groups.
  Web-based visualization:
    scp_output              Generate output files for single cell portal.
  MISC:
    check_indexes           Check CITE-Seq/hashing indexes to avoid index collision.

Options:
  -h, --help          Show help information.
  -v, --version       Show version.

Description:
  This is a tool for analyzing millions of single cell RNA-Seq data.

"""

from docopt import docopt
from docopt import DocoptExit
from . import __version__ as VERSION
from . import commands


def main():
    args = docopt(__doc__, version=VERSION, options_first=True)

    command_name = args["<command>"]
    command_args = args["<args>"]
    command_args.insert(0, command_name)

    try:
        command_class = getattr(commands, command_name)
    except AttributeError:
        print("Unknown command {cmd}!".format(cmd=command_name))
        raise DocoptExit()

    command = command_class(command_args)
    command.execute()


if __name__ == "__main__":
    main()
