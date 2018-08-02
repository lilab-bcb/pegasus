import os
import time
from .Base import Base
from ..annotate_cluster import run_annotate_cluster

class AnnotateCluster(Base):
    """
Annotate potential cell types for each cluster.

Usage:
  scrtools annotate_cluster [--json-file <file> --minimum-report-score <score> --do-not-use-non-de-genes] <input_h5ad_file> <output_file>
  scrtools annotate_cluster -h

Arguments:
  input_h5ad_file        Single cell data with DE analysis done by scrtools de_analysis.
  output_file            Output annotation file.

Options:
  --json-file <file>                    JSON file for markers. Could also be human_immune/mouse_immune/mouse_brain, which triggers scrtools to markers included in the package. [default: human_immune]
  --minimum-report-score <score>        Minimum cell type score to report a potential cell type. [default: 0.5]
  --do-not-use-non-de-genes             Do not count non DE genes as down-regulated.
  -h, --help                            Print out help information.

Outputs:
  output_file        This is a text file. For each cluster, all its putative cell types are listed in descending order of the cell type score. For each putative cell type, all markers support this cell type are listed. If one putative cell type has cell subtypes, all subtypes will be listed under this cell type.
  
Examples:
  scrtools annotate_cluster manton_bm_de.h5ad manton_bm.anno.txt
    """

    def execute(self):
        run_annotate_cluster(self.args['<input_h5ad_file>'], self.args['<output_file>'], float(self.args['--minimum-report-score']), ignoreNA = self.args['--do-not-use-non-de-genes'], json_file = self.args['--json-file'])
        time.sleep(1) # wait 1s to release the backed h5ad file
        