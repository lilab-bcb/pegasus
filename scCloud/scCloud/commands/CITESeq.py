from .Base import Base
from ..cite_seq import merge_rna_and_adt_data

class CITESeq(Base):
    """
Merge RNA and ADT matrices for CITE-Seq data.

Usage:
  scCloud merge_rna_adt <input_raw_gene_bc_matrices_h5.h5> <input_adt_csv_file> <antibody_control_csv> <output_10x.h5>
  scCloud merge_rna_adt -h

Arguments:
  input_raw_gene_bc_matrices_h5.h5        Input raw RNA expression matrix in 10x hdf5 format.
  input_adt_csv_file                      Input ADT (antibody tag) count matrix in CSV format.
  antibody_control_csv                    A CSV file containing the IgG control information for each antibody.
  output_10x.h5                           Merged output file in 10x hdf5 format.

Options:
  -h, --help        Print out help information.

Outputs:
  output_10x.h5        Output file in 10x hdf5 format. This file contains two groups --- one is for RNAs and the other is for ADTs.
  
Examples:
  scCloud merge_rna_adt example_raw_h5.h5 example_adt.csv antibody_control.csv example_merged_raw_10x.h5
    """

    def execute(self):
        run_cite_seq_pipeline(self.args['<input_raw_gene_bc_matrices_h5.h5>'], self.args['<input_adt_csv_file>'], self.args['<antibody_control_csv>'], self.args['<output_10x.h5>'])
