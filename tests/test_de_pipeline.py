import unittest
import os
import shutil
import scCloud.commands
from .test_util import assert_files_equal, assert_excel_equal, is_running_in_docker


class TestDePipeline(unittest.TestCase):

    def tearDown(self):
        os.remove('test_de.xlsx')
        os.remove('test_de.h5ad')

    def test_de_analysis(self):
        # de_analysis modifies h5ad file
        shutil.copy(os.path.join('tests', 'scCloud-test-data', 'output', 'test_cluster.h5ad'), 'test_de.h5ad')
        cmd = scCloud.commands.de_analysis(
            ['de_analysis', 'test_de.h5ad', 'test_de.xlsx', '--fisher',
             '--mwu', '--roc', '--labels', 'leiden_labels'])
        cmd.execute()
        # if is_running_in_docker():
        #     assert_files_equal(self, os.path.join('tests', 'output', 'test_de.xlsx'), 'test_de.xlsx')
        #     assert_files_equal(self, os.path.join('tests', 'output', 'test_de.h5ad'), 'test_de.h5ad')
        # else:
        assert_excel_equal(self, os.path.join('tests', 'scCloud-test-data', 'output', 'test_de.xlsx'), 'test_de.xlsx')


if __name__ == '__main__':
    unittest.main()
