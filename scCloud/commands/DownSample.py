from .Base import Base
from scCloud.tools import down_sample

class DownSample(Base):
    """
Downsample reads from a 10x molecule_info.h5 file.

Usage:
  scCloud down_sample [--random-state <seed>] <molecule_info.h5> <total_reads> <sampled_reads> <raw_feature_output_matrix_name>
  scCloud down_sample -h

Arguments:
  molecule_info.h5                    Input molecule info file (10x V3).
  total_reads                         Total number of sequenced reads.
  sampled_reads                       Target number of down-sampled reads.
  raw_feature_output_matrix_name      Output raw feature matrix name.

Options:
  --random-state <seed>        Seed for random number generator. [default: 0]
  -h, --help                   Print out help information.

Outputs:
  raw_feature_output_matrix.scCloud.h5        Output raw feature matrix.

Examples:
  scCloud down_sample --random-state 0 molecule_info.h5 400000000 100000000 output_matrix
    """
    def execute(self):
        down_sample(self.args['<molecule_info.h5>'], self.args['<raw_feature_output_matrix_name>'], int(self.args['<total_reads>']), int(self.args['<sampled_reads>']), random_state = int(self.args['--random-state']))
