from .Base import Base
from ..tools import run_conversion

class PARQUET(Base):
    """
Generate a PARQUET file for web-based visualization.

Usage:
  scCloud parquet [options] <input_h5ad_file> <output_name>
  scCloud parquet -h

Arguments:
  input_h5ad_file        Analyzed single cell data in h5ad format.
  output_name            Name prefix for the parquet file.

Options:
  -p <number>, --threads <number>        Number of threads used to generate the PARQUET file. [default: 1]
  -h, --help                             Print out help information.

Outputs:
  output_name.parquet        Generated PARQUET file that contains metadata and expression levels for every gene.

Examples:
  scCloud parquet manton_bm.h5ad manton_bm
    """

    def execute(self):
        run_conversion(self.args['<input_h5ad_file>'], self.args['<output_name>'], nthreads = int(self.args['--threads']))
