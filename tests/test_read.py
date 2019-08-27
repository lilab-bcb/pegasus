import unittest

import sccloud as scc
from .test_util import assert_adata_equal

import shutil
import os


class TestRead(unittest.TestCase):
    def tearDown(self):
        os.path.exists("test.h5ad") and os.remove("test.h5ad")
        os.path.exists("test_obsm_compound.h5ad") and os.remove(
            "test_obsm_compound.h5ad"
        )

    def test_mtx_v2(self):
        adata = scc.read_input(
            "tests/scCloud-test-data/input/hgmm_1k_filtered_gene_bc_matrices/hg19/matrix.mtx"
        )
        self.assertEqual(adata.shape[0], 504)

    def test_mtx_v3(self):
        adata = scc.read_input(
            "tests/scCloud-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/matrix.mtx.gz"
        )
        self.assertEqual(adata.shape[0], 1046)

    def test_mtx_v2_dir(self):
        adata = scc.read_input(
            "tests/scCloud-test-data/input/hgmm_1k_filtered_gene_bc_matrices/hg19/"
        )
        self.assertEqual(adata.shape[0], 504)

    def test_mtx_v3_dir(self):
        adata = scc.read_input(
            "tests/scCloud-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/"
        )
        self.assertEqual(adata.shape[0], 1046)

    def test_read_write_h5ad(self):
        adata = scc.read_input(
            "tests/scCloud-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/"
        )
        scc.write_output(adata, "test.h5ad")
        adata2 = scc.read_input("test.h5ad")
        assert_adata_equal(self, adata, adata2)

    def test_read_write_old_5ad(self):
        adata = scc.read_input(
            "tests/scCloud-test-data/input/test_obsm_compound.h5ad"
        )
        scc.write_output(adata, "test.h5ad")
        adata2 = scc.read_input("test.h5ad")
        assert_adata_equal(self, adata, adata2)

    def test_read_write_old_5ad_backed(self):
        shutil.copy(
            "tests/scCloud-test-data/input/test_obsm_compound.h5ad",
            "test_obsm_compound.h5ad",
        )
        adata = scc.read_input("test_obsm_compound.h5ad", h5ad_mode="r+")
        scc.write_output(adata, "test_obsm_compound.h5ad")
        adata2 = scc.read_input("test_obsm_compound.h5ad")
        assert_adata_equal(self, adata, adata2)

    def test_read_write_old_5ad_backed_whitelist(self):
        shutil.copy(
            "tests/scCloud-test-data/input/test_obsm_compound.h5ad",
            "test_obsm_compound.h5ad",
        )
        adata = scc.read_input("test_obsm_compound.h5ad", h5ad_mode="r+")
        scc.write_output(adata, "test_obsm_compound.h5ad", whitelist=['obs'])
        adata2 = scc.read_input("test_obsm_compound.h5ad")
        assert_adata_equal(self, adata, adata2)


if __name__ == "__main__":
    unittest.main()
