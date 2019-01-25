from .Base import Base
from ..tools import run_conversion

class FindMarkers(Base):
    """
Find markers using gradient boosting.

Usage:
  scCloud find_markers <input_h5ad_file> <label_attr> <n_jobs>
  scCloud find_markers -h

Arguments:
  input_h5ad_file        Analyzed single cell data in h5ad format.
  label_attr             Attribute for classification.
  n_jobs                 Number of threads.

Options:
  -h, --help                             Print out help information.

Examples:
  scCloud find_markers manton_bm.h5ad louvain_labels 20
    """

    def execute(self):
        run_benchmark(self.args['<input_h5ad_file>'], self.args['<label_attr>'], int(self.args['<n_jobs>']))
