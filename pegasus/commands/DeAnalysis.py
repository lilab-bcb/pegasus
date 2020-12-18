import os
from .Base import Base
from pegasus.tools import run_de_analysis


class DeAnalysis(Base):
    """
Perform DE analysis. Calculate Mann-Whitney U test and AUROC values by default.

Usage:
  pegasus de_analysis [options] (--labels <attr>) <input_data_file> <output_spreadsheet>
  pegasus de_analysis -h

Arguments:
  input_data_file        Single cell data with clustering calculated. DE results would be written back.
  output_spreadsheet     Output spreadsheet with DE results.

Options:
  --labels <attr>                  Use <attr> as cluster labels.
  --condition <attr>               Compute DE between conditions (one vs rest) in each cluster label if specified.
  --de-key <key>                   Store DE results into varm with key = <key>. [default: de_res]
  -p <threads>                     Use <threads> threads. [default: 1]
  --t                              Calculate Welch's t-test.
  --fisher                         Calculate Fisher's exact test.
  --temp-folder <temp_folder>      Joblib temporary folder for memmapping numpy arrays.
  --alpha <alpha>                  Control false discovery rate at <alpha>. [default: 0.05]
  --ndigits <ndigits>              Round non p-values and q-values to <ndigits> after decimal point in the excel. [default: 3]

  --quiet                          Do not show detailed intermediate outputs.
  -h, --help                       Print out help information.

Outputs:
  input_data_file        DE results would be written back to the 'varm' field with name set by --de-key <key>.
  output_spreadsheet     An excel spreadsheet containing DE results. Each cluster has two tabs in the spreadsheet. One is for up-regulated genes and the other is for down-regulated genes. If DE was performed on conditions within each cluster. Each cluster will have number of conditions tabs and each condition tab contains two spreadsheet: up for up-regulated genes and down for down-regulated genes.

Examples:
  pegasus de_analysis -p 26 --labels louvain_labels --t --fisher manton_bm.zarr.zip manton_bm_de.xlsx
    """

    def execute(self):
        run_de_analysis(
            self.args["<input_data_file>"],
            self.args["<output_spreadsheet>"],
            self.args["--labels"],
            condition=self.args["--condition"],
            de_key=self.args["--de-key"],
            n_jobs=int(self.args["-p"]),
            t=self.args["--t"],
            fisher=self.args["--fisher"],
            temp_folder=self.args["--temp-folder"],
            alpha=float(self.args["--alpha"]),
            ndigits=int(self.args["--ndigits"]),
            verbose=not self.args["--quiet"],
        )
