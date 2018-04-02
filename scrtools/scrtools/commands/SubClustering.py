from .Base import Base
from ..clustering import subcluster

class SubClustering(Base):
	"""
Run scanpy to obtain subclusters.

Usage:
  scrtools subcluster [options] <input_h5ad_file> <cluster_id> <output_name>
  scrtools subcluster -h

Arguments:
  input_h5ad_file        Single cell data with clustering done by Scanpy in h5ad file format.
  cluster_id             Tell us which cluster you want to perform subcluster task. ID starts from 1 and must be one of the cluster in input_h5ad_file.
  output_name            Output name. output_name.h5ad and output_name_var.h5ad will be generated.

Options:
  -p <number>, --threads <number>                  Number of threads. [default: 1]
  --correct-batch-effect                           Correct for batch effects for subclustering task. 
  
  --output-loom                                    Output loom-formatted file.
  --plot-by-side <attribute>                       Color tSNE by <attribute> and put it at the right side of louvain clusters.
  --legend-on-data                                 Put legends on data.

  --louvain-resolution <resolution>                Resolution parameter for the louvain clustering algorithm. [default: 1.3]
  --output-diagnosis-plots                         Output some diagnosis plots.

  -h, --help       Print out help information.

Examples:
  scrtools subcluster -p 20 --correct-batch-effect manton_bm.h5ad 1 manton_bm_1
	"""

	def execute(self):
		kwargs = {
			'nthreads' : int(self.args['--threads']),
			'correct_batch' : self.args['--correct-batch-effect'],
			'output_loom' : self.args['--output-loom'],
			'key' : self.args['--plot-by-side'],
			'legend_on_data' : self.args['--legend-on-data'],
			'resolution' : float(self.args['--louvain-resolution']),
			'diagnosis' : self.args['--output-diagnosis-plots']
		}

		subcluster(self.args['<input_h5ad_file>'], self.args['<cluster_id>'], self.args['<output_name>'], **kwargs)
