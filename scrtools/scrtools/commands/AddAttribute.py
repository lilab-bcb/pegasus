from .Base import Base
from ..preprocessing import add_attribute

class AddAttribute(Base):
	"""
Add one attribute to input.attr.csv file.

Usage:
  scrtools add_attribute <csv_file> <expression>
  scrtools add_attribute -h

Arguments:
  csv_file         scrtools aggregate_matrix generated attribute csv file. 
  expression       new_attribute:value or new_attribute:attribute=value or new_attribute:attr1+attr2+....

Options:
  -h, --help       Print out help information.

Examples:
  scrtools add_attribute manton_bm_cb.attr.csv batch:Individual=8
  scrtools add_attribute manton_bm_cb.attr.csv batch:Source+Individual
	"""

	def execute(self):
		add_attribute(self.args['<csv_file>'], self.args['<expression>'])
