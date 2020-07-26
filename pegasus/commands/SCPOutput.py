from .Base import Base
from pegasusio import read_input, write_output


class SCPOutput(Base):
    """
Generate outputs for single cell portal.

Usage:
  pegasus scp_output [--dense --round-to <ndigit>] <input_data_file> <output_name>
  pegasus scp_output -h

Arguments:
  input_data_file        Analyzed single cell data in zarr format.
  output_name            Name prefix for all outputted files.

Options:
  --dense                Output dense expression matrix instead.
  --round-to <ndigit>    Round expression to <ndigit> after the decimal point. [default: 2]
  -h, --help             Print out help information.

Outputs:
  output_name.scp.metadata.txt, output_name.scp.barcodes.tsv, output_name.scp.genes.tsv, output_name.scp.matrix.mtx, output_name.scp.*.coords.txt, output_name.scp.expr.txt         Files that single cell portal needs.

Examples:
  pegasus scp_output manton_bm.zarr.zip manton_bm
    """

    def execute(self):
        data = read_input(self.args["<input_data_file>"])
        write_output(
            data,
            self.args["<output_name>"],
            file_type = "scp",
            is_sparse = not self.args["--dense"],
            precision = int(self.args["--round-to"]),
        )
