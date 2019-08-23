import os
import unittest

import h5py

import sccloud as scc


class TestAggregate(unittest.TestCase):
    def tearDown(self):
        os.path.exists("aggregate_test.h5sc") and os.remove(
            "aggregate_test.h5sc"
        )

    def test_aggregate_10x_matrices(self):
        m1 = scc.read_input(
            "tests/scCloud-test-data/input/heart_1k_v3/filtered_feature_bc_matrix.h5",
            genome="mm10",
        )
        m2 = scc.read_input(
            "tests/scCloud-test-data/input/heart_1k_v2/filtered_gene_bc_matrices_h5.h5",
            genome="mm10",
        )
        scc.aggregate_matrices(
            "tests/scCloud-test-data/input/aggregate_test.csv",
            restrictions=[],
            attributes=["Version"],
            what_to_return="aggregate_test",
            google_cloud=False,
            select_singlets=False,
            ngene=None,
        )

        result = scc.read_input("aggregate_test.h5sc", genome="mm10")
        self.assertEqual(
            m1.shape[0] + m2.shape[0], result.shape[0], "Cell dimension is incorrect"
        )
        self.assertEqual(m1.shape[1], result.shape[1], "Gene dimension is incorrect")
        self.assertTrue(result.obs.get("Version") is not None, "Version not added")

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
        scc.aggregate_matrices(
            "tests/scCloud-test-data/input/aggregate_multi_genome.csv",
            restrictions=[],
            attributes=None,
            what_to_return="aggregate_test",
            google_cloud=False,
            select_singlets=False,
            ngene=None,
        )

        f = h5py.File("aggregate_test.h5sc", "r")
        self.assertIsNotNone(f["GRCh38"], "Genome not found")
        self.assertIsNotNone(f["mm10"], "Genome not found")
        with self.assertRaises(KeyError):
            f["mm9"]
        f.close()


if __name__ == "__main__":
    unittest.main()
