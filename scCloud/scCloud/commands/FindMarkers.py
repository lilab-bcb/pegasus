from .Base import Base
from scCloud.tools import run_find_markers

class FindMarkers(Base):
    """
Find markers using gradient boosting.

Usage:
  scCloud find_markers [options] <input_h5ad_file> <output_spreadsheet>
  scCloud find_markers -h

Arguments:
  input_h5ad_file        Single cell data after running the de_analysis.
  output_spreadsheet     Output spreadsheet with LightGBM detected markers.

Options:
  --labels <attr>              <attr> used as cluster labels. [default: louvain_labels]
  --remove-ribo                Remove ribosomal genes with either RPL or RPS as prefixes.
  --min-gain <gain>            Only report genes with a feature importance score (in gain) of at least <gain>. [default: 1.0]
  --random-state <seed>        Random state for initializing LightGBM and KMeans. [default: 0]
  -p <threads>                 Use <threads> threads. [default: 1]

  -h, --help                   Print out help information.

Outputs:
  output_spreadsheet     An excel spreadsheet containing detected markers. Each cluster has one tab in the spreadsheet and each tab has six columns, listing markers that are strongly up-regulated, weakly up-regulated, down-regulated and their associated LightGBM gains.

Examples:
  scCloud find_markers --labels louvain_labels --remove-ribo --min-gain 10.0 -p 10 manton_bm.h5ad manton_bm.markers.xlsx
    """

    def execute(self):
        run_find_markers(self.args['<input_h5ad_file>'], self.args['<output_spreadsheet>'], self.args['--labels'], n_jobs = int(self.args['-p']), \
            min_gain = float(self.args['--min-gain']), random_state = int(self.args['--random-state']), remove_ribo = self.args['--remove-ribo'])
