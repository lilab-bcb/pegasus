from .Base import Base
from ..annotate_cluster import run_annotate_cluster

class AnnotateCluster(Base):
	"""
Annotate potential cell types for each cluster.

Usage:
  scrtools annotate_cluster [--minimum-report-score <score>] <input_de_h5ad_file> <output_name>
  scrtools annotate_cluster -h

Arguments:
  input_de_h5ad_file        Single cell data with DE analysis done by scrtools de_analysis.
  output_name               Output name. output_name_fisher.anno.txt and output_name_t.anno.txt will be generated.

Options:
  --minimum-report-score <score>        Minimum cell type score to report a potential cell type. [default: 0.5]
  --do-not-use-non-de-genes             Do not count non DE genes as down regulated.
  -h, --help                            Print out help information.

Examples:
  scrtools annotate_cluster manton_bm_de.h5ad manton_bm
	"""

	def execute(self):
		run_annotate_cluster(self.args['<input_de_h5ad_file>'], self.args['<output_name>'], float(self.args['--minimum-report-score'], self.args['--do-not-use-non-de-genes']))
