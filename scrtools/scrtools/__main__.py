"""
Single cell RNA-Seq tools.

Usage:
  scrtools <command> [<args>...]
  scrtools -h | --help
  scrtools -v | --version

Sub-commands:
  Preprocessing related:
	aggregate_matrix        Aggregate hdf5-formatted 10x channel matrices into one big matrix. 
	add_attribute           Add one attribute for each channel.
	merge_matrix            Merge several 10x-hdf5-formatted matrices into one big matrix.

Options:
  -h, --help          Show help information.
  -v, --version       Show version.

Description:
  This is a tool for analyzing millions of single cell RNA-Seq data.

"""

from docopt import docopt
from docopt import DocoptExit
from . import __version__ as VERSION
from . import commands

def main():
	args = docopt(__doc__, version = VERSION, options_first = True)

	command_name = args['<command>']
	command_args = args['<args>']
	command_args.insert(0, command_name)

	try:
		command_class = getattr(commands, command_name)
	except AttributeError:
		print('Unknown command {cmd}!'.format(cmd = command_name))
		raise DocoptExit()
	
	command = command_class(command_args)
	command.execute()   




if __name__ == '__main__':
	main()
