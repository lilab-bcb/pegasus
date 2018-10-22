import os
import atexit

class Logging:
	def __init__(self, log_file, start_new = False):
		self.outputs = set()

		if start_new or not os.path.exists(log_file):
			self.fout = open(log_file, 'w')
		else:
			with open(log_file) as fin:
				for line in fin:
					self.outputs.add(line.strip())
			self.fout = open(log_file, 'a')

		atexit.register(self.close)

	def add_output(self, output_file):
		if output_file not in self.outputs:
			self.outputs.add(output_file)
			self.fout.write(output_file + '\n')

	def close(self):
		self.fout.close()
