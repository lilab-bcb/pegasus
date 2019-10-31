from .Base import Base
from pegasus.tools import aggregate_matrices


class AggregateMatrix(Base):
    """
Aggregate 10x matrices from each channel into one big matrix.

Usage:
  pegasus aggregate_matrix <csv_file> <output_name> [--restriction <restriction>... --attributes <attributes> --genome <genome> --select-only-singlets --minimum-number-of-genes <ngene>]
  pegasus aggregate_matrix -h

Arguments:
  csv_file          Input csv-formatted file containing information of each scRNA-Seq run. Each row must contain at least 2 columns - Sample, sample name and Location, location of the channel-specific count matrix in either 10x v2/v3, DGE, mtx, csv or loom format. If matrix is in DGE, mtx or csv format, an additional Reference column is required. Additionally an optional RenamedReference column can be used to rename a reference.
  output_name       The output file name.

Options:
  --restriction <restriction>...           Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value or name:~value,..,value, where ~ refers to not. You can specifiy multiple restrictions by setting this option multiple times.
  --attributes <attributes>                Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.
  --genome <genome>                        Default genome to use.
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
            genome=self.args["--genome"],
            select_singlets=self.args["--select-only-singlets"],
            ngene=self.convert_to_int(self.args["--minimum-number-of-genes"]),
        )
