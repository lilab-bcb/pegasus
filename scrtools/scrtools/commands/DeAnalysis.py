from .Base import Base
from ..de_analysis import run_de_analysis

class DeAnalysis(Base):
	"""
Perform t-test and Fisher's exact test.

Usage:
  scrtools de_analysis [--fold-change <threshold>] <input_h5ad_file> <output_name>
  scrtools de_analysis -h

Arguments:
  input_h5ad_file        Single cell data with clustering done by Scanpy in h5ad file format.
  output_name            Output name. output_name_de.h5ad, output_name_de_analysis_fisher.xlsx, and output_name_de_analysis_t.xlsx will be generated.

Options:
  --fold-change <threshold>        Minimum fold change in either percentage (fisher test) or log expression (t test) to report a DE gene. [default: 1.5]
  -h, --help                       Print out help information.

Examples:
  scrtools de_analysis manton_bm.h5ad manton_bm
	"""

	def execute(self):
		run_de_analysis(self.args['<input_h5ad_file>'], self.args['<output_name>'], float(self.args['--fold-change']))
