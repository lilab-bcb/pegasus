import os
import unittest

import h5py

import scCloud as sc


class TestAggregate(unittest.TestCase):

    def test_aggregate_10x_matrices(self):
        m1 = sc.tools.read_input('data/heart_1k_v3/filtered_gene_bc_matrices_h5.h5', genome='mm10',
                                 mode='r')
        m2 = sc.tools.read_input('data/heart_1k_v2/filtered_gene_bc_matrices_h5.h5', genome='mm10',
                                 mode='r')
        sc.tools.aggregate_10x_matrices('data/aggregate_test.csv', restrictions=[], attributes=['Version'],
                                        output_file='test.h5',
                                        google_cloud=False, select_singlets=False, ngene=None, is_dropseq=None)

        result = sc.tools.read_input('test.h5', genome='mm10', mode='r')
        self.assertEqual(m1.shape[0] + m2.shape[0], result.shape[0], "Cell dimension is incorrect")
        self.assertEqual(m1.shape[1], result.shape[1], "Gene dimension is incorrect")
        self.assertTrue(result.obs.get('Version') is not None, 'Version not added')

        m1_result = result[list(range(m1.shape[0])), :]
        m2_result = result[list(range(m1.shape[0], m1.shape[0] + m2.shape[0])), :]
        self.assertEqual((m1_result.X != m1.X).sum(), 0, 'Values differ')
        self.assertEqual((m2_result.X != m2.X).sum(), 0, 'Values differ')
        self.assertTrue(m1_result.obs.index.values[0].startswith('heart_1k_v3'), 'Prefix not added')
        self.assertTrue(m2_result.obs.index.values[0].startswith('heart_1k_v2'), 'Prefix not added')

    def tearDown(self):
        os.remove('test.h5')

    def test_multi_genome(self):
        sc.tools.aggregate_10x_matrices('data/aggregate_multi_genome.csv', restrictions=[], attributes=None,
                                        output_file='test.h5', google_cloud=False, select_singlets=False,
                                        ngene=None, is_dropseq=None)

        f = h5py.File('test.h5', 'r')
        self.assertIsNotNone(f['GRCh38'], 'Genome not found')
        self.assertIsNotNone(f['mm10'], 'Genome not found')
        with self.assertRaises(KeyError):
            f['mm9']
        f.close()
