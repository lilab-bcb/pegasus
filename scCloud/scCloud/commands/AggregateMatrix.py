from .Base import Base
from ..tools import aggregate_10x_matrices, Logging

class AggregateMatrix(Base):
	"""
Aggregate 10x matrices from each channel into one big matrix.

Usage:
  scCloud aggregate_matrix <csv_file> <output_name> [--genome <genome> --restriction <restriction>... --attributes <attributes> --google-cloud --input-type <type>]
  scCloud aggregate_matrix -h

Arguments:
  csv_file          Input csv-formatted file containing information of each 10x channel. Each row must contain at least 3 columns --- Sample, sample name; Location, location of the channel-specific count matrix in 10x format (e.g. /sample/filtered_gene_bc_matrices_h5.h5); Reference, genome reference used for 10x cellranger.
  output_name       The output file name. 

Options:
  --genome <genome>                Genome reference. [default: GRCh38]
  --restriction <restriction>...   Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value or name:~value,..,value, where ~ refers to not. You can specifiy multiple restrictions by setting this option multiple times.
  --attributes <attributes>        Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.
  --google-cloud                   If files are stored in google cloud. Assuming google cloud sdk is installed.
  --input-type <type>              Input type, could be either 'gene', 'dropseq' or 'ADT'. [default: gene]
  -h, --help                       Print out help information.

Outputs:
  output_name_10x.h5        A 10x-formatted HDF5 file containing the count matrix and associated attributes.
       
Examples:
  scCloud aggregate_matrix --genome GRCh38 --restriction Source:BM,CB --restriction Individual:1-8 --attributes Source,Platform Manton_count_matrix.csv manton_bm_cb
	"""

	def execute(self):
		aggregate_10x_matrices(self.args['<csv_file>'], 
			self.args['--genome'], 
			self.args['--restriction'], 
			self.split_string(self.args['--attributes']), 
			self.args['<output_name>'],
			google_cloud = self.args['--google-cloud'],
      input_type = self.args['--input-type'])
