from .Base import Base
# from ..clustering import cluster, cluster2

class Clustering(Base):
	"""
Run scanpy to obtain top level clusters.

Usage:
  scrtools cluster [options] <input_name> <output_name>
  scrtools cluster -h

Arguments:
  input_name       Input name of 10x format gene count matrix. Please make sure that input_name_10x.h5 and input_name.attr.csv exist.
  output_name      Output name. output_name.h5ad and output_name_var.h5ad will be generated.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --genome <genome>                                Genome name. [default: GRCh38]
  
  --output-filtration-results <spreadsheet>        Output filtration results into <spreadsheet>.
  --output-loom                                    Output loom-formatted file.
  --correct-batch-effect                           Correct for batch effects.
  --batch-group-by <groupby>                       Group batches by. If this option is off, make all batches as a single group.
  
  --import-attributes <attributes>                 Import attributes contained in the comma-separated list into the analysis object.
  --plot-by-side <attribute>                       Color tSNE by <attribute> and put it at the right side of louvain clusters.
  --legend-on-data                                 Put legends on data.

  --min-genes <number>                             Only keep cells with at least <number> of genes. [default: 500]
  --max-genes <number>                             Only keep cells with less than <number> of genes. [default: 6000]
  --mito-prefix <prefix>                           Prefix for mitochondrial genes. [default: MT-]
  --percent-mito <ratio>                           Only keep cells with mitochondrial ratio less than <ratio>. [default: 0.1]

  --gene-percent-cells <ratio>                     Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes. [default: 0.0005]

  --counts-per-cell-after <number>                 Total counts per cell after normalization. [default: 1e5]
  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --output-diagnosis-plots                         Output some diagnosis plots.

  --input-h5ad                                     Input is one h5ad file.
  --nPC <number>                                   Number of PCs.

  -h, --help       Print out help information.

Examples:
  scrtools cluster -p 20 --correct-batch-effect manton_bm manton_bm
	"""

	def execute(self):
		kwargs = {
			'nthreads' : int(self.args['--threads']),
			'genome' : self.args['--genome'],
			'filtr_xlsx' : self.args['--output-filtration-results'],
			'output_loom' : self.args['--output-loom'],
			'correct_batch' : self.args['--correct-batch-effect'],
			'groupby' : self.split_string(self.args['--batch-group-by']),
			'attrs' : self.split_string(self.args['--import-attributes']),
			'key' : self.args['--plot-by-side'],
			'legend_on_data' : self.args['--legend-on-data'],
			'min_genes' : int(self.args['--min-genes']),
			'max_genes' : int(self.args['--max-genes']),
			'mito_prefix' : self.args['--mito-prefix'],
			'percent_mito' : float(self.args['--percent-mito']),
			'percent_cells' : float(self.args['--gene-percent-cells']),
			'norm_count' : float(self.args['--counts-per-cell-after']),
			'resolution' : float(self.args['--louvain-resolution']),
			'diagnosis' : self.args['--output-diagnosis-plots'],
			'input_h5ad' : self.args['--input-h5ad'],
			'nPC' : self.args['--nPC']
		}

		# if kwargs['nPC'] is not None:
		# 	kwargs['nPC'] = int(kwargs['nPC'])
		# 	cluster2(self.args['<input_name>'], self.args['<output_name>'], **kwargs)
		# else:
		# 	cluster(self.args['<input_name>'], self.args['<output_name>'], **kwargs)
