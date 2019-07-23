import os
from .Base import Base
from ..annotate_cluster import run_annotate_cluster, annotate_anndata_object

class AnnotateCluster(Base):
	"""
Annotate potential cell types for each cluster. This command has two forms: the first form generates putative annotations and the second form write annotations into the h5ad object.

Usage:
  scCloud annotate_cluster [--json-file <file> --minimum-report-score <score> --do-not-use-non-de-genes] <input_h5ad_file> <output_file>
  scCloud annotate_cluster --annotation <annotation_string> <input_h5ad_file>
  scCloud annotate_cluster -h

Arguments:
  input_h5ad_file        Single cell data with DE analysis done by scCloud de_analysis.
  output_file            Output annotation file.

Options:
  --json-file <file>                      JSON file for markers. Could also be human_immune/mouse_immune/mouse_brain/human_brain, which triggers scCloud to markers included in the package. [default: human_immune]
  --minimum-report-score <score>          Minimum cell type score to report a potential cell type. [default: 0.5]
  --do-not-use-non-de-genes               Do not count non DE genes as down-regulated.

  --annotation <annotation_string>        Write cell type annotations in <annotation_string> into <input_h5ad_file>. <annotation_string> has this format: 'anno_attr:anno_1;anno_2;...;anno_n'. 'anno_attr' is the annotation attribute in the h5ad object and anno_i is the annotation for cluster i.

  -h, --help                              Print out help information.

Outputs:
  output_file        This is a text file. For each cluster, all its putative cell types are listed in descending order of the cell type score. For each putative cell type, all markers support this cell type are listed. If one putative cell type has cell subtypes, all subtypes will be listed under this cell type.
  
Examples:
  scCloud annotate_cluster manton_bm.h5ad manton_bm.anno.txt
  scCloud annotate_cluster --annotation "anno:T cells;B cells;NK cells;Monocytes" manton_bm.h5ad
    """

	def execute(self):
		if self.args['<output_file>'] is not None:
			run_annotate_cluster(self.args['<input_h5ad_file>'], self.args['<output_file>'], float(self.args['--minimum-report-score']), ignoreNA = self.args['--do-not-use-non-de-genes'], json_file = self.args['--json-file'])
		else:
			annotate_anndata_object(self.args['<input_h5ad_file>'], self.args['--annotation'])
