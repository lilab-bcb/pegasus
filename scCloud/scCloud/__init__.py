__version__ = "0.10.0"

try:
	get_ipython
except NameError:
	import matplotlib
	matplotlib.use('Agg')

from . import tools, plotting, demuxEM, cite_seq, misc
