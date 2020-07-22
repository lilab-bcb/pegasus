from .Base import Base
from pegasus.check_sample_indexes import run_check_sample_indexes


class CheckSampleIndexes(Base):
    """
Check for index collision between 10x scRNA-seq index sets and CITE-Seq/hashing indexes. 

This command can also be used to find the maximum number of mismatches allowed among HTO/ADT barcodes.

Usage:
  pegasus check_indexes [--num-mismatch <mismatch> --report <n_report>] <index_file>
  pegasus check_indexes -h

Arguments:
  index_file        Index file containing CITE-Seq/hashing index sequences. One sequence per line. Multiple columns are allowed. But columns must be separated by comma and the first column must be the index sequence.

Options:
  --num-mismatch <mismatch>        Number of mismatch allowed for each index sequence. [default: 1]
  --report <n_report>              Number of valid 10x GA indexes to report. Default is not to calculate valid GA indexes. [default: -1]
  -h, --help                       Print out help information.

Outputs:
  If --report is not set, <index_file> should include all scRNA-seq/CITE-seq/hashing indexes. This program first report the minimum hamming distance between any pair of indexes and also the maximum number of mismatches can be set [(hamming-dist - 1) // 2]. If the maximum number of mismatches is smaller than <mismatch>, an index collision error message will be generated.

  If --report is set, assume <index_file> only contain CITE-seq/hashing indexes. If there is no index collision within <index_file>, up to <n_report> number of valid 10x scRNA-seq indexes will be printed to the standard output.

Examples:
  pegasus check_indexes --num-report 8 index_file.txt
    """

    def execute(self):
        run_check_sample_indexes(
            self.args["<index_file>"],
            n_mis=int(self.args["--num-mismatch"]),
            n_report=int(self.args["--report"]),
        )
