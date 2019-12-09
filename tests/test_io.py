import os
import shutil
import unittest

import numpy as np
import pandas as pd

import pegasus as pg
from .test_util import assert_adata_equal


class TestRead(unittest.TestCase):
    def tearDown(self):
        for f in ["test.csv", "test.h5ad", "test_obsm_compound.h5ad"]:
            os.path.exists(f) and os.remove(f)
        os.path.exists('test_mtx') and shutil.rmtree('test_mtx')

    def test_mtx_v2(self):
        adata = pg.read_input(
            "tests/pegasus-test-data/input/hgmm_1k_filtered_gene_bc_matrices/hg19/matrix.mtx"
        )
        self.assertEqual(adata.shape[0], 504)

    def test_mtx_v3(self):
        adata = pg.read_input(
            "tests/pegasus-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/matrix.mtx.gz"
        )
        self.assertEqual(adata.shape[0], 1046)

    def test_mtx_v2_dir(self):
        adata = pg.read_input(
            "tests/pegasus-test-data/input/hgmm_1k_filtered_gene_bc_matrices/hg19/"
        )
        self.assertEqual(adata.shape[0], 504)

    def test_write_mtx(self):
        adata = pg.read_input(
            "tests/pegasus-test-data/input/heart_1k_v3/filtered_feature_bc_matrix.h5"
        )
        adata.var['test'] = 1.0
        adata.obs['test'] = 1.0
        output_dir = 'test_mtx/mm10'
        pg.write_output(adata, os.path.join(output_dir, 'matrix.mtx.gz'))
        adata2 = pg.read_input(output_dir)
        del adata2.obs['Channel']  # get channel from csv
        adata2.obs = adata2.obs.join(pd.read_csv(os.path.join(output_dir, 'obs.csv.gz'), index_col=0))
        adata2.var = adata2.var.join(pd.read_csv(os.path.join(output_dir, 'var.csv.gz'), index_col=0))
        del adata2.var['featuretype']
        assert_adata_equal(self, adata, adata2, obs_blacklist=['Channel'])

    def test_mtx_v3_dir(self):
        adata = pg.read_input(
            "tests/pegasus-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/"
        )
        self.assertEqual(adata.shape[0], 1046)

    def test_csv(self):
        df = pd.DataFrame(index=["a", "b", "c"], data=dict(a=[1, 2, 3], b=[4, 5, 6]))
        df.to_csv("test.csv")
        adata = pg.read_input("test.csv", genome="test").T
        np.testing.assert_array_equal(df.values, adata.X.toarray())
        np.testing.assert_array_equal(df.index.values, adata.obs.index.values)
        np.testing.assert_array_equal(df.columns.values, adata.var.index.values)

        for chunk_size in [1, 2, 3, 4]:
            adata_chunks = pg.read_input(
                "test.csv", genome="test", chunk_size=chunk_size
            ).T
            assert_adata_equal(self, adata, adata_chunks)

    def test_read_write_h5ad(self):
        adata = pg.read_input(
            "tests/pegasus-test-data/input/hgmm_1k_v3_filtered_feature_bc_matrix/"
        )
        pg.write_output(adata, "test.h5ad")
        adata2 = pg.read_input("test.h5ad")
        assert_adata_equal(self, adata, adata2)

    def test_read_write_old_5ad(self):
        adata = pg.read_input("tests/pegasus-test-data/input/test_obsm_compound.h5ad")
        pg.write_output(adata, "test.h5ad")
        adata2 = pg.read_input("test.h5ad")
        assert_adata_equal(self, adata, adata2)

    def test_read_write_old_5ad_backed(self):
        shutil.copy(
            "tests/pegasus-test-data/input/test_obsm_compound.h5ad",
            "test_obsm_compound.h5ad",
        )
        adata = pg.read_input("test_obsm_compound.h5ad", h5ad_mode="r+")
        pg.write_output(adata, "test_obsm_compound.h5ad")
        adata2 = pg.read_input("test_obsm_compound.h5ad")
        assert_adata_equal(self, adata, adata2)

    def test_read_write_old_5ad_backed_whitelist(self):
        shutil.copy(
            "tests/pegasus-test-data/input/test_obsm_compound.h5ad",
            "test_obsm_compound.h5ad",
        )
        adata = pg.read_input("test_obsm_compound.h5ad", h5ad_mode="r+")
        pg.write_output(adata, "test_obsm_compound.h5ad", whitelist=["obs"])
        adata2 = pg.read_input("test_obsm_compound.h5ad")
        assert_adata_equal(self, adata, adata2)


if __name__ == "__main__":
    unittest.main()
