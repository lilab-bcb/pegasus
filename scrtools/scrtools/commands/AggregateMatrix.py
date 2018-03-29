from .Base import Base
from ..preprocessing import aggregate_10x_matrices

class AggregateMatrix(Base):
	"""
Aggregate 10x matrices from each channel into one big matrix.

Usage:
  scrtools aggregate_matrix <csv_file> <output_name> [--genome <genome> --restriction <restriction>... --attributes <attributes>]
  scrtools aggregate_matrix -h

Arguments:
  csv_file          Input csv-formatted file containing information of each 10x channel. Each row must contain at least 3 columns --- Sample, sample name; Location, folder of 10x Run; Reference, genome reference used for 10x cellranger. scrtools will look for {Location}/{Sample}/filtered_gene_bc_matrices_h5.h5.
  output_name       The output file prefix. Two files will be generated: the aggregated data matrix, output_name_10x.h5, and the channel attribute file, output_name.attr.csv.

Options:
  --genome <genome>                Genome reference. [default: GRCh38]
  --restriction <restriction>...   Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. You can specifiy multiple restrictions by setting this option multiple times.
  --attributes <attributes>        Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.
  -h, --help                       Print out help information.

Examples:
  scrtools aggregate_matrix --genome GRCh38 --restrictions Source:BM,CB;Individual:1-8 --attributes Source,Platform Manton_count_matrix.csv manton_bm_cb
	"""

	def parse_attributes(self, attributes):
		return attributes.split(',') if attributes is not None else []

	def execute(self):
		aggregate_10x_matrices(self.args['<csv_file>'], 
			self.args['--genome'], 
			self.args['--restriction'], 
			self.split_string(self.args['--attributes']), 
			self.args['<output_name>'])
