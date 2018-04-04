from .Base import Base
from ..plotting import make_plots

class Plotting(Base):
	"""
Generate cluster composition plots.

Usage:
  scrtools plot composition [options] --attribute <attr> <input_h5ad_file> <output_file>
  scrtools plot diffmap [--dim <dim>] --attributes <attrs> <input_h5ad_file> <output_file>
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

  --attributes <attrs>      A comma-separated list of attributes to plot side-by-side for diffmap.
  --dim <dim>               Diffusion map dimension, can be either '2d' or '3d'. [default: 2d]
  -h, --help                Print out help information.

Examples:
  scrtools plot --attribute Individual manton_bm.h5ad manton_bm.png
	"""

	def execute(self):
		kwargs = {
			'attr' : self.args['--attribute'],
			'style' : self.args['--style'],
			'stacked' : not self.args['--not-stacked'],
			'logy' : self.args['--log-y'],
			'sizes' : [float(x) for x in self.args['--sizes'].split(',')],
			'rmove' : float(self.args['--rmove']),
			'wshrink' : float(self.args['--shrink']),

			'attrs' : self.split_string(self.args['--attributes']),
			'dim' : self.args['--dim']
		}

		keyword = ''
		if self.args['composition']:
			keyword = 'composition'
		elif self.args['diffmap']:
			keyword = 'diffmap'
		else:
			assert False

		make_plots(self.args['<input_h5ad_file>'], keyword, self.args['<output_file>'], **kwargs)
