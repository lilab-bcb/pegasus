from .Base import Base
from pegasus.tools import aggregate_matrices


class AggregateMatrix(Base):
    """
Aggregate sample count matrices from each channel into one big matrix.

Usage:
  pegasus aggregate_matrix <csv_file> <output_name> [--restriction <restriction>... --attributes <attributes> --default-reference <reference> --select-only-singlets --minimum-number-of-genes <ngene>]
  pegasus aggregate_matrix -h

Arguments:
  csv_file          Input csv-formatted file containing information of each sc/snRNA-seq sample. This file must contain at least 2 columns - Sample, sample name and Location, location of the sample count matrix in either 10x v2/v3, DGE, mtx, csv, tsv or loom format. Additionally, an optional Reference column can be used to select samples generated from a same reference (e.g. mm10). If the count matrix is in either DGE, mtx, csv, tsv, or loom format, the value in this column will be used as the reference since the count matrix file does not contain reference name information. In addition, the Reference column can be used to aggregate count matrices generated from different genome versions or gene annotations together under a unified reference. For example, if we have one matrix generated from mm9 and the other one generated from mm10, we can write mm9_10 for these two matrices in their Reference column. Pegasus will change their references to 'mm9_10' and use the union of gene symbols from the two matrices as the gene symbols of the aggregated matrix. For HDF5 files (e.g. 10x v2/v3), the reference name contained in the file does not need to match the value in this column. In fact, we use this column to rename references in HDF5 files. For example, if we have two HDF files, one generated from mm9 and the other generated from mm10. We can set these two files' Reference column value to 'mm9_10', which will rename their reference names into mm9_10 and the aggregated matrix will contain all genes from either mm9 or mm10. This renaming feature does not work if one HDF5 file contain multiple references (e.g. mm10 and GRCh38).
  output_name       The output file name.

Options:
  --restriction <restriction>...           Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value or name:~value,..,value, where ~ refers to not. You can specifiy multiple restrictions by setting this option multiple times.
  --attributes <attributes>                Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.
  --default-reference <reference>          If sample count matrix is in either DGE, mtx, csv, tsv or loom format and there is no Reference column in the csv_file, use <reference> as the reference.
  --select-only-singlets                   If we have demultiplexed data, turning on this option will make pegasus only include barcodes that are predicted as singlets.
  --minimum-number-of-genes <ngene>        Only keep barcodes with at least <ngene> expressed genes.
  -h, --help                               Print out help information.

Outputs:
  output_name.h5sc        A pegasus-formatted HDF5 file containing the count matrices and associated attributes.

Examples:
  pegasus aggregate_matrix --restriction Source:BM,CB --restriction Individual:1-8 --attributes Source,Platform Manton_count_matrix.csv manton_bm_cb
    """

    def execute(self):
        aggregate_matrices(
            self.args["<csv_file>"],
            what_to_return=self.args["<output_name>"],
            restrictions=self.args["--restriction"],
            attributes=self.split_string(self.args["--attributes"]),
            default_ref=self.args["--default-reference"],
            select_singlets=self.args["--select-only-singlets"],
            ngene=self.convert_to_int(self.args["--minimum-number-of-genes"]),
        )
