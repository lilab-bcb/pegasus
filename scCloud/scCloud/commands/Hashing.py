from .Base import Base
from ..pipeline import run_hashing_pipeline

class Hashing(Base):
    """
Run the hashing pipeline for cell-hashing/nuclei-hashing.

Usage:
  scCloud hashing --hash-type <type> [options] <input_adt_csv_file> <input_raw_gene_bc_matrices_h5.h5> <output_name>
  scCloud hashing -h

Arguments:
  input_adt_csv_file                      Input ADT (antibody tag) count matrix in CSV format.
  input_raw_gene_bc_matrices_h5.h5        Input raw RNA expression matrix in 10x hdf5 format.
  output_name                             Output name. All outputs will use it as the prefix.

Options:
  --hash-type <type>                         The hash type of the data. <type> can be 'cell' for cell-hashing and 'nuclei' for nuclei-hashing.

  -p <number>, --threads <number>            Number of threads. [default: 1]
  --genome <genome>                          Reference genome name. If not provided, we will infer it from the expression matrix file.

  --min-num-genes <number>                   We only demultiplex cells/nuclei with at least <number> expressed genes. [default: 100]
  --max-background-probability <prob>        Any cell/nucleus with no less than <prob> background probability will be marked as unknown. [default: 0.8]
  --prior-on-samples <prior>                 The sparse prior put on samples. [default: 0.0]

  -h, --help                                 Print out help information.

Outputs:
  output_name.h5ad                Demultiplexing results in h5ad format.
  output_name.demux_10x.h5        RNA expression matrix with demultiplexed sample identities in 10x's hdf5 format.
  
Examples:
  scCloud hashing -p 8 --hash-type cell sample_adt.csv sample_raw_gene_bc_matrices_h5.h5 sample_demux
    """

    def execute(self):
        kwargs = {
            'hash_type' : self.args['--hash-type'],
            'n_jobs' : int(self.args['--threads']),
            'genome' : self.args['--genome'],
            'min_num_genes' : int(self.args['--min-num-genes']),
            'unknown_rate' : float(self.args['--max-background-probability']),
            'alpha_value' : float(self.args['--prior-on-samples'])
        }

        run_cite_seq_pipeline(self.args['<input_adt_csv_file>'], self.args['<input_raw_gene_bc_matrices_h5.h5>'], self.args['<output_name>'], **kwargs)
