from .Base import Base
from ..tools import run_scp_output

class SCPOutput(Base):
    """
Generate outputs for single cell portal.

Usage:
  scCloud scp_output [--dense --round-to <ndigit>] <input_h5ad_file> <output_name>
  scCloud scp_output -h

Arguments:
  input_h5ad_file        Analyzed single cell data in h5ad format.
  output_name            Name prefix for all outputted files.

Options:
  --dense                Output dense expression matrix instead.
  --round-to <ndigit>    Round expression to <ndigit> after the decimal point. [default: 2]
  -h, --help             Print out help information.

Outputs:
  output_name.scp.metadata.txt, output_name.scp.barcodes.tsv, output_name.scp.genes.tsv, output_name.scp.matrix.mtx, output_name.scp.*.coords.txt, output_name.scp.expr.txt         Files that single cell portal needs.

Examples:
  scCloud scp_output manton_bm.h5ad manton_bm
    """

    def execute(self):
        run_scp_output(self.args['<input_h5ad_file>'], self.args['<output_name>'], not self.args['--dense'], int(self.args['--round-to']))
