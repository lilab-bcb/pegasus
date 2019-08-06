import unittest

import scCloud as sc
from .test_util import assert_adata_equal


class TestRead(unittest.TestCase):

    def test_mtx_v2(self):
        adata = sc.io.read_input('tests/scCloud-test-data/input/hgmm_1k_filtered_gene_bc_matrices/hg19/matrix.mtx')
        self.assertEqual(adata.shape[0], 504)

    def test_mtx_v3(self):
        adata = sc.io.read_input('tests/scCloud-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/matrix.mtx.gz')
        self.assertEqual(adata.shape[0], 1046)

    def test_mtx_v2_dir(self):
        adata = sc.io.read_input('tests/scCloud-test-data/input/hgmm_1k_filtered_gene_bc_matrices/hg19/')
        self.assertEqual(adata.shape[0], 504)

    def test_mtx_v3_dir(self):
        adata = sc.io.read_input('tests/scCloud-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/')
        self.assertEqual(adata.shape[0], 1046)

    def test_read_write_h5ad(self):
        adata = sc.io.read_input('tests/scCloud-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/')
        sc.io.write_output(adata, 'test.h5ad')
        adata2 = sc.io.read_input('test.h5ad', h5ad_mode='a')
        assert_adata_equal(self, adata, adata2)


if __name__ == '__main__':
    unittest.main()
