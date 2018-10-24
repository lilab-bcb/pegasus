from .Base import Base
from ..pipeline import run_demuxEM_pipeline

class DemuxEM(Base):
    """
Run the demuxEM pipeline for cell-hashing/nuclei-hashing data.

Usage:
  scCloud demuxEM --hash-type <type> [options] <input_adt_csv_file> <input_raw_gene_bc_matrices_h5.h5> <output_name>
  scCloud demuxEM -h

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
  --prior-on-samples <prior>                 The sparse prior put on samples.
  --random-state <seed>                      The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]

  --generate-diagnostic-plots                Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc.
  --generate-gender-plot <genes>             Generate violin plots using gender-specific genes (e.g. Xist). <gene> is a comma-separated list of gene names. 

  -h, --help                                 Print out help information.

Outputs:
  output_name.demux_10x.h5                               RNA expression matrix with demultiplexed sample identities in 10x's hdf5 format.
  output_name_ADTs.h5ad                                  Antibody tag matrix in h5ad format.
  output_name_demux.h5ad                                 Demultiplexed RNA count matrix in h5ad format.
  output_name.ambient_hashtag.hist.png                   Optional output. A histogram plot depicting hashtag distributions of empty droplets and non-empty droplets.
  output_name.background_probabilities.bar.png           Optional output. A bar plot visualizing the estimated hashtag background probability distribution.
  output_name.real_content.hist.png                      Optional output. A histogram plot depicting hashtag distributions of not-real-cells and real-cells as defined by total number of expressed genes in the RNA assay.
  output_name.rna_demux.hist.png                         Optional output. A histogram plot depicting RNA UMI distribution for singlets, doublets and unknown cells.
  output_name.gene_name.violin.png                       Optional outputs. Violin plots depicting gender-specific gene expression across samples. We can have multiple plots if a gene list is provided in '--generate-gender-plot' option.
  
Examples:
  scCloud demuxEM -p 8 --hash-type cell --generate-diagnostic-plots sample_adt.csv sample_raw_gene_bc_matrices_h5.h5 sample_output
    """

    def execute(self):
        kwargs = {
            'hash_type' : self.args['--hash-type'],
            'n_jobs' : int(self.args['--threads']),
            'genome' : self.args['--genome'],
            'min_num_genes' : int(self.args['--min-num-genes']),
            'unknown_rate' : float(self.args['--max-background-probability']),
            'alpha_value' : self.args['--prior-on-samples'],
            'random_state' : int(self.args['--random-state']),
            'gen_plots' : self.args['--generate-diagnostic-plots'],
            'gen_gender_plot' : self.split_string(self.args['--generate-gender-plot'])
        }

        run_demuxEM_pipeline(self.args['<input_adt_csv_file>'], self.args['<input_raw_gene_bc_matrices_h5.h5>'], self.args['<output_name>'], **kwargs)
