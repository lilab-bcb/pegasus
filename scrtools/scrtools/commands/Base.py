from docopt import docopt
from docopt import DocoptExit

class Base:
	"""Base class for the commands"""

	def __init__(self, command_args):
		self.args = docopt(self.__doc__, argv = command_args)

	def split_string(self, astring, sep = ','):
		return astring.split(sep) if astring is not None else None

	def execute(self):
		raise NotImplementedError
