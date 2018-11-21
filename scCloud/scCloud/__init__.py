__version__ = "0.8.0"

try:
	get_ipython
except NameError:
	import matplotlib
	matplotlib.use('Agg')

from . import tools
from . import plotting
from . import misc
from . import demuxEM
