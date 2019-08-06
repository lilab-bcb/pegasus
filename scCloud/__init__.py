__version__ = "0.13.0"

try:
    get_ipython
except NameError:
    import matplotlib

    matplotlib.use("Agg")

from . import tools, plotting, demuxEM, cite_seq, misc

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
