import os
import unittest

import anndata
import h5py
import numpy as np
import pandas as pd
import pegasus as pg


class TestAggregate(unittest.TestCase):
    def tearDown(self):
        for file in ["aggregate_test.h5sc", "test1_dge.txt.gz", "test2_dge.txt.gz", "aggregate_test"]:
            os.path.exists(file) and os.remove(file)

    def test_aggregate_union(self):
        m1 = anndata.AnnData(X=np.zeros((1, 3)), var=pd.DataFrame(index=['a', 'b', 'c']),
            obs=pd.DataFrame(index=['c1']))
        pd.DataFrame(columns=m1.obs.index, index=m1.var.index, data=m1.X.T).to_csv('test1_dge.txt.gz', index_label='id',
            sep='\t')
        m2 = anndata.AnnData(X=np.zeros((1, 3)), var=pd.DataFrame(index=['b', 'd', 'e']),
            obs=pd.DataFrame(index=['c1']))
        pd.DataFrame(columns=m2.obs.index, index=m2.var.index, data=m2.X.T).to_csv('test2_dge.txt.gz', index_label='id',
            sep='\t')
        pd.DataFrame(
            data=dict(Sample=['s1', 's2'], Location=['test1_dge.txt.gz', 'test2_dge.txt.gz'],
                Reference=['a', 'a'])).to_csv(
            'aggregate_test.csv', index=False)
        result = pg.aggregate_matrices(
            "aggregate_test.csv",
            what_to_return="AnnData"
        )
        self.assertEqual(result.shape[0], 2)
        self.assertEqual(result.shape[1], 5)

    def test_aggregate_10x_matrices(self):
        m1 = pg.read_input(
            "tests/pegasus-test-data/input/heart_1k_v3/filtered_feature_bc_matrix.h5",
            genome="mm10",
        )
        m2 = pg.read_input(
            "tests/pegasus-test-data/input/heart_1k_v2/filtered_gene_bc_matrices_h5.h5",
            genome="mm10",
        )
        pg.aggregate_matrices(
            "tests/pegasus-test-data/input/aggregate_test.csv",
            what_to_return='aggregate_test',
        )

        result = pg.read_input("aggregate_test.h5sc", genome="mm10")
        self.assertEqual(
            m1.shape[0] + m2.shape[0], result.shape[0], "Cell dimension is incorrect"
        )
        self.assertEqual(m1.shape[1], result.shape[1], "Feature dimension is incorrect")

        m1_result = result[list(range(m1.shape[0])), :]
        m2_result = result[list(range(m1.shape[0], m1.shape[0] + m2.shape[0])), :]
        self.assertEqual((m1_result.X != m1.X).sum(), 0, "Values differ")
        self.assertEqual((m2_result.X != m2.X).sum(), 0, "Values differ")
        self.assertTrue(
            m1_result.obs.index.values[0].startswith("heart_1k_v3"), "Prefix not added"
        )
        self.assertTrue(
            m2_result.obs.index.values[0].startswith("heart_1k_v2"), "Prefix not added"
        )

    def test_multi_genome(self):
        pg.aggregate_matrices(
            "tests/pegasus-test-data/input/aggregate_multi_genome.csv",
            what_to_return='aggregate_test',
        )

        f = h5py.File("aggregate_test.h5sc", "r")
        self.assertIsNotNone(f["GRCh38"], "Genome not found")
        self.assertIsNotNone(f["mm10"], "Genome not found")
        with self.assertRaises(KeyError):
            f["mm9"]
        f.close()


if __name__ == "__main__":
    unittest.main()
