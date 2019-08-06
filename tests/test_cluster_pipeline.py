import unittest
import os
import scCloud.commands
from .test_util import assert_adata_files_equal


class TestClusterPipeline(unittest.TestCase):
    def tearDown(self):
        os.remove('test_cluster.h5ad')
        os.remove('test_cluster.hvg.pdf')

    def test_cluster(self):
        cmd = scCloud.commands.cluster(
            ['cluster', os.path.join('tests', 'scCloud-test-data', 'input', '3k_pbmc'),
             'test_cluster', '--run-leiden',
             '--run-approximated-leiden', '--run-tsne', '--run-umap',
             '--run-net-tsne', '--run-net-fitsne', '--run-net-umap', '--run-fitsne', '--run-fle',
             '--run-net-fle', '--plot-hvg', '--run-louvain', '--run-approximated-louvain'])
        cmd.execute()
        #
        # if is_running_in_docker:
        #     assert_files_equal(self, os.path.join('tests', 'output', 'test_cluster.h5ad'), 'test_cluster.h5ad')
        #     assert_files_equal(self, os.path.join('tests', 'output', 'test_cluster.hvg.pdf'), 'test_cluster.hvg.pdf')
        # else:
        # TODO diff pdfs
        assert_adata_files_equal(self, os.path.join('tests', 'scCloud-test-data', 'output', 'test_cluster.h5ad'),
            'test_cluster.h5ad')


if __name__ == '__main__':
    unittest.main()
