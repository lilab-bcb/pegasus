from .Base import Base
from pegasus.tools import run_find_markers


class FindMarkers(Base):
    """
Find markers using gradient boosting.

Usage:
  pegasus find_markers [options] <input_data_file> <output_spreadsheet>
  pegasus find_markers -h

Arguments:
  input_data_file        Single cell data after running the de_analysis.
  output_spreadsheet     Output spreadsheet with LightGBM detected markers.

Options:
  -p <threads>                 Use <threads> threads. [default: 1]
  --labels <attr>              <attr> used as cluster labels. [default: louvain_labels]
  --de-key <key>               Key for storing DE results in 'varm' field. [default: de_res]
  --remove-ribo                Remove ribosomal genes with either RPL or RPS as prefixes.
  --min-gain <gain>            Only report genes with a feature importance score (in gain) of at least <gain>. [default: 1.0]
  --random-state <seed>        Random state for initializing LightGBM and KMeans. [default: 0]

  -h, --help                   Print out help information.

Outputs:
  output_spreadsheet     An excel spreadsheet containing detected markers. Each cluster has one tab in the spreadsheet and each tab has six columns, listing markers that are strongly up-regulated, weakly up-regulated, down-regulated and their associated LightGBM gains.

Examples:
  pegasus find_markers --labels louvain_labels --remove-ribo --min-gain 10.0 -p 10 manton_bm.zarr.zip manton_bm.markers.xlsx
    """

    def execute(self):
        run_find_markers(
            self.args["<input_data_file>"],
            self.args["<output_spreadsheet>"],
            self.args["--labels"],
            de_key=self.args["--de-key"],
            n_jobs=int(self.args["-p"]),
            min_gain=float(self.args["--min-gain"]),
            random_state=int(self.args["--random-state"]),
            remove_ribo=self.args["--remove-ribo"],
        )
