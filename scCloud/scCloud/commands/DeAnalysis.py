import os
from .Base import Base
from ..tools import run_de_analysis

class DeAnalysis(Base):
    """
Perform DE analysis.

Usage:
  scCloud de_analysis [--labels <attr> --subset <attr:value> -p <threads> --alpha <alpha> --fisher --mwu --roc --temp-folder <temp_folder>] <input_h5ad_file> <output_spreadsheet>
  scCloud de_analysis -h

Arguments:
  input_h5ad_file        Single cell data with clustering calculated. DE results would be written back.
  output_spreadsheet     Output spreadsheet with DE results.

Options:
  --labels <attr>                  <attr> used as cluster labels. [default: louvain_labels]
  --subset <attr:value>            Perform DE analysis only on a subset of cells with attr == value. 
  --alpha <alpha>                  Control false discovery rate at <alpha>. [default: 0.05]
  --fisher                         Calculate Fisher's exact test.
  --mwu                            Calculate Mann-Whitney U test.
  --roc                            Calculate area under cuver in ROC curve.
  --temp-folder <temp_folder>      Joblib temporary folder for memmapping numpy arrays.
  -p <threads>                     Use <threads> threads. [default: 1]
  
  -h, --help                       Print out help information.

Outputs:
  input_h5ad_file        DE results would be written back to the 'var' fields, provided --subset option is not set.
  output_spreadsheet     An excel spreadsheet containing DE results. Each cluster has two tabs in the spreadsheet. One is for up-regulated genes and the other is for down-regulated genes.
  output_h5ad_file       Only present if --subset option is set. The file name is output_spreadsheet - [.xlsx] + [.h5ad]

Examples:
  scCloud de_analysis --labels louvain_labels -p 26 --fisher --mwu --roc manton_bm.h5ad manton_bm_de.xlsx
    """

    def execute(self):
        run_de_analysis(self.args['<input_h5ad_file>'], self.args['<output_spreadsheet>'], self.args['--labels'], int(self.args['-p']), float(self.args['--alpha']), self.args['--fisher'], self.args['--mwu'], self.args['--roc'], self.args['--subset'], self.args['--temp-folder'])
