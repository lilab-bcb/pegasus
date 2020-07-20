import os
from .Base import Base
from pegasus.tools import run_de_analysis


class DeAnalysis(Base):
    """
Perform DE analysis.

Usage:
  pegasus de_analysis [options] <input_data_file> <output_spreadsheet>
  pegasus de_analysis -h

Arguments:
  input_data_file        Single cell data with clustering calculated. DE results would be written back.
  output_spreadsheet     Output spreadsheet with DE results.

Options:
  -p <threads>                     Use <threads> threads. [default: 1]
  --labels <attr>                  <attr> used as cluster labels. [default: louvain_labels]
  --result-key <key>               Store DE results into AnnData varm with key = <key>. [default: de_res]
  --auc                            Calculate area under ROC (AUROC) and area under Precision-Recall (AUPR).
  --t                              Calculate Welch's t-test.
  --fisher                         Calculate Fisher's exact test.
  --mwu                            Calculate Mann-Whitney U test.
  --temp-folder <temp_folder>      Joblib temporary folder for memmapping numpy arrays.
  --alpha <alpha>                  Control false discovery rate at <alpha>. [default: 0.05]
  --ndigits <ndigits>              Round non p-values and q-values to <ndigits> after decimal point in the excel. [default: 3]

  --quiet                          Do not show detailed intermediate outputs.
  -h, --help                       Print out help information.

Outputs:
  input_data_file        DE results would be written back to the 'varm' field with name set by --result-key <key>.
  output_spreadsheet     An excel spreadsheet containing DE results. Each cluster has two tabs in the spreadsheet. One is for up-regulated genes and the other is for down-regulated genes.

Examples:
  pegasus de_analysis -p 26 --labels louvain_labels --auc --t --fisher --mwu manton_bm.zarr.zip manton_bm_de.xlsx
    """

    def execute(self):
        run_de_analysis(
            self.args["<input_data_file>"],
            self.args["<output_spreadsheet>"],
            self.args["--labels"],
            result_key=self.args["--result-key"],
            n_jobs=int(self.args["-p"]),
            auc=self.args["--auc"],
            t=self.args["--t"],
            fisher=self.args["--fisher"],
            mwu=self.args["--mwu"],
            temp_folder=self.args["--temp-folder"],
            verbose=not self.args["--quiet"],
            alpha=float(self.args["--alpha"]),
            ndigits=int(self.args["--ndigits"]),
        )
