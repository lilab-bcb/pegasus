from .Base import Base
from ..check_sample_indexes import run_check_sample_indexes

class CheckSampleIndexes(Base):
    """
Check for index collision between 10x' index sets and CITE-Seq/hashing indexes.

Usage:
  scCloud check_indexes [--num-mismatch <mismatch> --num-report <report>] <index_file>
  scCloud check_indexes -h

Arguments:
  index_file        Index file containing CITE-Seq/hashing index sequences. One sequence per line.

Options:
  --num-mismatch <mismatch>        Number of mismatch allowed for each index sequence. [default: 1]
  --num-report <report>            Number of valid 10x indexes to report. Default is to report all valid indexes. [default: 9999]
  -h, --help                       Print out help information.

Outputs:
  Up to <report> number of valid 10x indexes will be printed out to standard output. 

Examples:
  scCloud check_indexes --num-report 8 index_file.txt
    """

    def execute(self):
        run_check_sample_indexes(self.args['<index_file>'], n_mis = int(self.args['--num-mismatch']), n_report = int(self.args['--num-report']))
