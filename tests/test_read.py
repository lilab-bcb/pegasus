import unittest

import scCloud as sc


class TestRead(unittest.TestCase):

	def test_mtx_v2(self):
		adata = sc.tools.read_input('data/hgmm_1k_filtered_gene_bc_matrices/hg19/matrix.mtx')
		self.assertEqual(adata.shape[0], 504)

	def test_mtx_v3(self):
		adata = sc.tools.read_input('data/hgmm_1k_v3_filtered_feature_bc_matrix/matrix.mtx.gz')
		self.assertEqual(adata.shape[0], 1046)

	def test_mtx_v2_dir(self):
		adata = sc.tools.read_input('data/hgmm_1k_filtered_gene_bc_matrices/hg19/')
		self.assertEqual(adata.shape[0], 504)

	def test_mtx_v3_dir(self):
		adata = sc.tools.read_input('data/hgmm_1k_v3_filtered_feature_bc_matrix/')
		self.assertEqual(adata.shape[0], 1046)
