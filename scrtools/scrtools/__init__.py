__version__ = "0.3.0"

try:
	get_ipython
except NameError:
	import matplotlib
	matplotlib.use('Agg')

from . import tools
from . import plotting
from . import misc
