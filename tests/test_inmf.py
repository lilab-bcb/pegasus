import unittest

import numpy as np

class TestINMF(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestINMF, self).__init__(*args, **kwargs)
        self.n_cells = 1043
        self.n_features = 18952
        self.n_hvfs = 2000
        self.n_factors = 20
        self.n_batches = 2

    def test_h5ad(self):
        from anndata import read_h5ad
        adata = read_h5ad("tests/inmf_result.mm10-rna.h5ad")

        self.assertEqual(adata.shape, (self.n_cells, self.n_features), "Count matrix shape not correct!")
        self.assertEqual(adata.uns['H'].shape, (self.n_cells, self.n_factors), "H shape not correct!")
        self.assertEqual(adata.uns['V'].shape, (self.n_batches, self.n_factors, self.n_hvfs), "V shape not correct!")
        self.assertEqual(adata.uns['W'].shape, (self.n_hvfs, self.n_factors), "W shape not correct!")
        self.assertEqual(adata.obsm['X_inmf'].shape, (self.n_cells, self.n_factors), "iNMF embedding shape not correct!")

    def test_loom(self):
        import loompy
        with loompy.connect("tests/inmf_result.mm10-rna.loom") as ds:
            self.assertEqual(ds.shape, (self.n_features, self.n_cells), "Count matrix shape not correct!")
            self.assertEqual(ds.ca['X_inmf'].shape, (self.n_cells, self.n_factors), "iNMF embedding shape not correct!")

    def test_zarr(self):
        import pegasusio as io
        data = io.read_input("tests/inmf_result.zarr.zip")

        self.assertEqual(data.shape, (self.n_cells, self.n_features), "Count matrix shape not correct!")
        self.assertEqual(data.uns['H'].shape, (self.n_cells, self.n_factors), "H shape not correct!")
        self.assertEqual(data.uns['V'].shape, (self.n_batches, self.n_factors, self.n_hvfs), "V shape not correct!")
        self.assertEqual(data.uns['W'].shape, (self.n_hvfs, self.n_factors), "W shape not correct!")

        self.assertEqual(data.obsm['X_inmf'].shape, (self.n_cells, self.n_factors), "iNMF embedding shape not correct!")


if __name__ == "__main__":
    unittest.main()
