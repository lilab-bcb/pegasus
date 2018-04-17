from .Base import Base
from ..plotting import make_interactive_plots

class iPlotting(Base):
	"""
Generate cluster composition plots.

Usage:
  scrtools iplot diffmap [options] --attribute <attr> <input_h5ad_file> <output_html_file>
  scrtools iplot -h

Arguments:
  input_h5ad_file        Single cell data with clustering done by Scanpy in h5ad file format.
  output_html_file       Output interactive htl plot file name.

Options:
  --attribute <attr>        Use attribute <attr> as labels in the plot.
  --real                    <attr> is real valued.
  --log10                   If take log10 of real values.

  -h, --help             Print out help information.

Examples:
  scrtools iplot diffmap --attribute annotation manton_bm.h5ad manton_bm_diffmap.html
	"""
	def execute(self):
		kwargs = {
			'attr' : self.args['--attribute'],
			'real' : self.args['--real'],
			'log10' : self.args['--log10']
		}

		keyword = ''
		if self.args['diffmap']:
			keyword = 'diffmap'
		else:
			assert False

		make_interactive_plots(self.args['<input_h5ad_file>'], keyword, self.args['<output_html_file>'], **kwargs)
