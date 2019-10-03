import os
from .Base import Base
from pegasus.misc import show_attributes


class View(Base):
    """
Show sample and gene attributes contained in the h5ad file.

Usage:
  pegasus view [--show-attributes --show-gene-attributes --show-values-for-attributes <attributes>] <input_h5ad_file>
  pegasus view -h

Arguments:
  input_h5ad_file        Analyzed single cell data in h5ad format.

Options:
  --show-attributes                                Show the available sample attributes in the input dataset.
  --show-gene-attributes                           Show the available gene attributes in the input dataset.
  --show-values-for-attributes <attributes>        Show the available values for specified attributes in the input dataset. <attributes> should be a comma-separated list of attributes.
  -h, --help                                       Print out help information.

Examples:
  pegasus view --show-attributes manton_bm.h5ad
  pegasus view --show-gene-attributes manton_bm.h5ad
  pegasus view --show-values-for-attributes louvain_labels,Condition manton_bm.h5ad
    """

    def execute(self):
        show_attributes(
            self.args["<input_h5ad_file>"],
            self.args["--show-attributes"],
            self.args["--show-gene-attributes"],
            self.args["--show-values-for-attributes"],
        )
