import os
from .Base import Base
from pegasus.annotate_cluster import run_annotate_cluster, annotate_data_object


class AnnotateCluster(Base):
    """
Annotate potential cell types for each cluster. This command has two forms: the first form generates putative annotations and the second form write annotations into the data object.

Usage:
  pegasus annotate_cluster [--marker-file <file> --de-test <test> --de-alpha <alpha> --de-key <key> --minimum-report-score <score> --do-not-use-non-de-genes] <input_data_file> <output_file>
  pegasus annotate_cluster --annotation <annotation_string> <input_data_file>
  pegasus annotate_cluster -h

Arguments:
  input_data_file        Single cell data with DE analysis done by pegasus de_analysis.
  output_file            Output annotation file.

Options:
  --marker-file <file>                    JSON file for markers. Could also be human_immune/mouse_immune/human_brain/mouse_brain/human_lung, which triggers pegasus to markers included in the package. [default: human_immune]
  --de-test <test>                        DE test to use to infer cell types. [default: t]
  --de-alpha <alpha>                      False discovery rate to control family-wise error rate. [default: 0.05]
  --de-key <key>                          Keyword where the DE results store in varm. [default: de_res]
  --minimum-report-score <score>          Minimum cell type score to report a potential cell type. [default: 0.5]
  --do-not-use-non-de-genes               Do not count non DE genes as down-regulated.

  --annotation <annotation_string>        Write cell type annotations in <annotation_string> into <input_h5ad_file>. <annotation_string> has this format: 'anno_name:clust_name:anno_1;anno_2;...;anno_n'. 'anno_name' is the annotation attribute in the h5ad object, 'clust_name' is the attribute with cluster ids, and anno_i is the annotation for cluster i.

  -h, --help                              Print out help information.

Outputs:
  output_file        This is a text file. For each cluster, all its putative cell types are listed in descending order of the cell type score. For each putative cell type, all markers support this cell type are listed. If one putative cell type has cell subtypes, all subtypes will be listed under this cell type.

Examples:
  pegasus annotate_cluster manton_bm.zarr.zip manton_bm.anno.txt
  pegasus annotate_cluster --annotation "anno:louvain_labels:T cells;B cells;NK cells;Monocytes" manton_bm.zarr.zip
    """

    def execute(self):
        if self.args["<output_file>"] is not None:
            run_annotate_cluster(
                self.args["<input_data_file>"],
                self.args["<output_file>"],
                self.args["--marker-file"],
                de_test=self.args["--de-test"],
                de_alpha=float(self.args["--de-alpha"]),
                de_key=self.args["--de-key"],
                threshold=float(self.args["--minimum-report-score"]),
                ignore_nonde=self.args["--do-not-use-non-de-genes"],
            )
        else:
            annotate_data_object(
                self.args["<input_data_file>"], self.args["--annotation"]
            )
