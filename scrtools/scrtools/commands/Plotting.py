from .Base import Base
from ..plotting import make_plots

class Plotting(Base):
	"""
Generate cluster composition plots.

Usage:
  scrtools plot [options] --attribute <attr> <input_h5ad_file> <output_file>
  scrtools plot -h

Arguments:
  input_h5ad_file        Single cell data with clustering done by Scanpy in h5ad file format.
  output_file            Output image file with suffix from .pdf, .png.

Options:
  --attribute <attr>        Plot <attr> against cluster labels.
  --style <style>           Plot styles. Can be either 'frequency', 'count', or 'normalized'. [default: frequency]
  --not-stacked             Do not stack bars.
  --log-y                   Plot y axis in log10 scale.
  --sizes <sizes>           Figure sizes in inches, w x h, separated by comma. [default: 6,4]
  --rmove <rmove>           Move legend to right for <rmove>. [default: 1.1]
  --shrink <shrink>         Shrink figures. [default: 0.8]
  -h, --help                Print out help information.

Examples:
  scrtools plot --attribute Individual manton_bm.h5ad manton_bm.png
	"""

	def execute(self):
		kwargs = {
			'attr' : self.args['--attribute'],
			'style' : self.args['--style'],
			'stack' : not self.args['--not-stacked'],
			'logy' : self.args['--log-y']
			'sizes' : [float(x) for x in self.args['--sizes'].split(',')],
			'rmove' : float(self.args['--rmove']),
			'wshrink' : float(self.args['--shrink'])
		}
		print(kwargs)
		return
		make_plots(self.args['<input_h5ad_file>'], 'composition', self.args['output_file'], kwargs)
