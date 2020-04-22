import pandas as pd

from .Base import Base


class Enrichment(Base):
    """
Perform enrichment analysis using gprofiler.

Usage:
  pegasus enrichment [options] <markers_spreadsheet> <output_spreadsheet>
  pegasus enrichment -h

Arguments:
  markers_spreadsheet        Output spreadsheet from de_analysis command.
  output_spreadsheet         File containing enrichment results.

Options:
  --organism <value>               Organism. See https://biit.cs.ut.ee/gprofiler/page/organism-list for full list. [default: hsapiens]
  --enrichment_threshold <value>   Include enrichment results with corrected p-value less than this threshold. [default: 0.05]
  --max_genes <value>              Maximum number of genes to include in query. [default: 100]
  -h, --help                       Print out help information.

Outputs:
  output     An xlsx file containing enrichment results.

Examples:
  pegasus enrichment manton_bm.de.xlsx manton_bm_enrichment.xlsx
    """

    def execute(self):
        d = pd.read_excel(self.args["<markers_spreadsheet>"], sheet_name=None)
        output_spreadsheet = self.args['<output_spreadsheet>']
        organism = self.args["--organism"]
        enrichment_threshold = float(self.args["--enrichment_threshold"])
        max_genes = int(self.args['--max_genes'])
        from gprofiler import GProfiler
        gp = GProfiler(return_dataframe=True)
        query = {}
        for key in d.keys():
            features = d[key]['feature'].values.tolist()
            query[key] = features[0:max_genes]

        result = gp.profile(organism=organism, query=query, user_threshold=enrichment_threshold)
        result.to_excel(output_spreadsheet, index=False)
