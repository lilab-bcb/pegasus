from pegasus.tools import run_conversion

from .Base import Base


class PARQUET(Base):
    """
Generate a PARQUET file for web-based visualization.

Usage:
  pegasus parquet [options] <input_h5ad_file> <output_name>
  pegasus parquet -h

Arguments:
  input_h5ad_file        Analyzed single cell data in h5ad format.
  output_name            Name prefix for the parquet file.

Options:
  -p <number>, --threads <number>           Number of threads used to generate the PARQUET file. [default: 1]
  -r <number>, --row_groups <number>        Number of parquet file row groups.
  -b, --backed                              Load h5ad file in backed mode instead of fully loading it into memory.
  -h, --help                                Print out help information.

Outputs:
  output_name.parquet        Generated PARQUET file that contains metadata and expression levels for every gene.

Examples:
  pegasus parquet -p 8 manton_bm.h5ad manton_bm
    """

    def execute(self):
        run_conversion(
            self.args["<input_h5ad_file>"],
            self.args["<output_name>"],
            nthreads=int(self.args["--threads"]),
            backed=self.args['--backed'],
            n_row_groups=int(self.args["--row_groups"]) if self.args.get("--row_groups", None) is not None else None
        )
