from .Base import Base
from pegasus.io import read_input, write_output


class ConvertH5SC(Base):
    """
Convert H5SC file into 10x V2 formatted HDF5 file.

Usage:
  pegasus convert_h5sc <input.h5sc> <output_10x.h5>
  pegasus convert_h5sc -h

Arguments:
  input.h5sc           Input file name.
  output_10x.h5        Output file name.

Options:
  -h, --help           Print out help information.

Outputs:
  output_name.h5sc        A pegasus-formatted HDF5 file containing the count matrices and associated attributes.

Examples:
  pegasus convert_h5sc input.h5sc output_10x.h5
    """

    def execute(self):
      data = read_input(self.args["<input.h5sc>"], return_type = "MemData")
      write_output(data, self.args["<output_10x.h5>"])

