from .Base import Base
from ..tools import aggregate_10x_matrices

class AggregateMatrix(Base):
	"""
Aggregate 10x matrices from each channel into one big matrix.

Usage:
  scCloud aggregate_matrix <csv_file> <output_name> [--restriction <restriction>... --attributes <attributes> --google-cloud --select-only-singlets --minimum-number-of-genes <ngene> --dropseq-genome <genome>]
  scCloud aggregate_matrix -h

Arguments:
  csv_file          Input csv-formatted file containing information of each scRNA-Seq run. Each row must contain at least 2 columns --- Sample, sample name and Location, location of the channel-specific count matrix in either 10x format (e.g. /sample/filtered_gene_bc_matrices_h5.h5) or dropseq format (e.g. /sample/sample.umi.dge.txt.gz).
  output_name       The output file name. 

Options:
  --restriction <restriction>...           Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value or name:~value,..,value, where ~ refers to not. You can specifiy multiple restrictions by setting this option multiple times.
  --attributes <attributes>                Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file.
  --google-cloud                           If files are stored in google cloud. Assuming google cloud sdk is installed.
  --select-only-singlets                   If we have demultiplexed data, turning on this option will make scCloud only include barcodes that are predicted as singlets.
  --minimum-number-of-genes <ngene>        Only keep barcodes with at least <ngene> expressed genes.
  --dropseq-genome <genome>                If inputs are dropseq data, this option needs to turn on and provides the reference genome name.

  -h, --help                               Print out help information.

Outputs:
  output_name_10x.h5        A 10x-formatted HDF5 file containing the count matrices and associated attributes.
       
Examples:
  scCloud aggregate_matrix --restriction Source:BM,CB --restriction Individual:1-8 --attributes Source,Platform Manton_count_matrix.csv manton_bm_cb
	"""

	def execute(self):
		aggregate_10x_matrices(self.args['<csv_file>'], 
			self.args['--restriction'], 
			self.split_string(self.args['--attributes']), 
			self.args['<output_name>'] + '_10x.h5',
			google_cloud = self.args['--google-cloud'],
      select_singlets = self.args['--select-only-singlets'],
      ngene = int(self.args['--minimum-number-of-genes']) if self.args['--minimum-number-of-genes'] is not None else None,
      is_dropseq = self.args['--dropseq-genome'])
