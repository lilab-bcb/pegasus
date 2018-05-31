import os
from .Base import Base
from ..annotate_cluster import run_annotate_cluster
from ..tools import Logging

class AnnotateCluster(Base):
    """
Annotate potential cell types for each cluster.

Usage:
  scrtools annotate_cluster [--minimum-report-score <score> --do-not-use-non-de-genes] <input_h5ad_file> <output_file>
  scrtools annotate_cluster -h

Arguments:
  input_h5ad_file        Single cell data with DE analysis done by scrtools de_analysis.
  output_file            Output annotation file.

Options:
  --minimum-report-score <score>        Minimum cell type score to report a potential cell type. [default: 0.5]
  --do-not-use-non-de-genes             Do not count non DE genes as down regulated.
  -h, --help                            Print out help information.

Examples:
  scrtools annotate_cluster manton_bm_de.h5ad manton_bm.anno.txt
    """

    def execute(self):
        run_annotate_cluster(self.args['<input_h5ad_file>'], self.args['<output_file>'], float(self.args['--minimum-report-score']), self.args['--do-not-use-non-de-genes'])

        logger = Logging(os.path.splitext(self.args['<input_h5ad_file>'])[0] + ".log")
        logger.add_output(self.args['<output_file>'])
