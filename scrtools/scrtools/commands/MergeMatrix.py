from .Base import Base
from ..preprocessing import merge_10x_matrices

class MergeMatrix(Base):
	"""
Merge several matrices into one big matrix.

Usage:
  scrtools merge_matrix --input <input_name>... <output_name> [--genome <genome> --symbols <symbols> --attributes <attributes>]
  scrtools merge_matrix -h

Arguments:
  output_name        Output matrix name. output_name_10x.h5 and output_name.attr.csv will be generated.

Options:
  -i <input_name>..., --input <input_name>...        Names of input matrices. input_name_10x.h5 and input_name.attr.csv must exist.
  --genome <genome>                                  Name of the genome. [default: GRCh38]
  --symbols <symbols>                                A comma-separated list of symbols representing each input matrix.
  --attributes <attributes>                          A comma-separated list of attributes. When merging matrices, the matrix symbol will be added in front of the attributes in the list.
  -h, --help                                         Print out help information.

Examples:
  scrtools merge_matrix --symbols BM,CB --attributes Individual -i manton_bm -i manton_cb manton_bm_cb
	"""

  

	def execute(self):
		merge_10x_matrices(self.args['--input'], 
			self.split_string(self.args['--symbols']),
			self.args['--genome'], 
			self.split_string(self.args['--attributes']), 
			self.args['<output_name>'])
