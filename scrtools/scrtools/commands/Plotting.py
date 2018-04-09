from .Base import Base
from ..plotting import make_static_plots

class Plotting(Base):
	"""
Generate cluster composition plots.

Usage:
  scrtools plot composition [options] --attribute <attr> <input_h5ad_file> <output_name>
  scrtools plot diffmap [options] --attributes <attr> <input_h5ad_file> <output_name>
  scrtools plot -h

Arguments:
  input_h5ad_file        Single cell data with clustering done by Scanpy in h5ad file format.
  output_name            Output image file name.

Options:
  --ext <ext>               Output file extension. [default: png]
  --attribute <attr>        Plot <attr> against cluster labels.

  --style <style>           Plot styles. Can be either 'frequency', 'count', or 'normalized'. [default: frequency]
  --not-stacked             Do not stack bars.
  --log-y                   Plot y axis in log10 scale.
  --sizes <sizes>           Figure sizes in inches, w x h, separated by comma. [default: 6,4]
  --rmove <rmove>           Move legend to right for <rmove>. [default: 1.1]
  --shrink <shrink>         Shrink figures. [default: 0.8]

  --attributes <attrs>                A comma-separated list of attributes to plot side-by-side for diffmap.
  --dim <dim>                         Diffusion map dimension, can be either '2d' or '3d'. [default: 2d]
  --angles <angles>                   Which angles you want to plot, could be 0, 1, 2, or all. For 2d plots, 0 = 1,2; 1 = 2,3; 2 = 1,3. For 3d plots, 0 = 1,2,3; 1 = 2,3,1; 2 = 1,3,2. [default: all]
  --point-size <size>                 Point size in the plots.
  --legend-fontsize <fontsize>.       Legend font size.

  -h, --help                Print out help information.

Examples:
  scrtools plot composition --attribute Individual manton_bm.h5ad manton_bm
	"""
	def execute(self):
		kwargs = {
			'ext' : self.args['--ext'],

			'attr' : self.args['--attribute'],
			'style' : self.args['--style'],
			'stacked' : not self.args['--not-stacked'],
			'logy' : self.args['--log-y'],
			'sizes' : [float(x) for x in self.args['--sizes'].split(',')],
			'rmove' : float(self.args['--rmove']),
			'wshrink' : float(self.args['--shrink']),

			'attrs' : self.split_string(self.args['--attributes']),
			'dim' : self.args['--dim'],
			'angles' : self.args['--angles'],
			'ptsize' : float(self.args['--point-size']) if self.args['--point-size'] is not None else None,
			'lfsize' : float(self.args['--legend-fontsize']) if self.args['--legend-fontsize'] is not None else None
		}

		keyword = ''
		if self.args['composition']:
			keyword = 'composition'
		elif self.args['diffmap']:
			keyword = 'diffmap'
		else:
			assert False

		make_static_plots(self.args['<input_h5ad_file>'], keyword, self.args['<output_name>'], **kwargs)
