import unittest
import subprocess
import os
import scCloud as sc
import pandas as pd
import numpy as np
import scipy.sparse


class TestPipeline(unittest.TestCase):
    def tearDown(self):
        os.remove('test_3k_pbmc.h5ad')
        os.remove('test_3k_pbmc.hvg.pdf')


    def test_cluster(self):
        subprocess.call(
            ['scCloud', 'cluster', os.path.join('data', '3k_pbmc'), 'test_3k_pbmc', '--run-leiden',
             '--run-approximated-leiden', '--run-tsne', '--run-umap',
             '--run-net-tsne', '--run-net-fitsne', '--run-net-umap', '--run-fitsne'])
        test_data = sc.tools.read_input('test_3k_pbmc.h5ad')
        data = sc.tools.read_input(os.path.join('output', 'test_3k_pbmc.h5ad'))
        self.assertEqual((test_data.X[()] != data.X[()]).sum(), 0)
        pd.testing.assert_frame_equal(test_data.obs, data.obs)
        pd.testing.assert_frame_equal(test_data.var, data.var)
        self.assertListEqual(list(test_data.uns.keys()), list(data.uns.keys()))
        self.assertListEqual(list(test_data.obsm_keys()), list(data.obsm_keys()))
        for key in data.uns_keys():
            test_val = test_data.uns[key]
            val = data.uns[key]
            if scipy.sparse.issparse(val):
                val = val.toarray()
                test_val = test_val.toarray()
            np.testing.assert_array_equal(test_val, val)
        for key in data.obsm_keys():
            test_val = test_data.obsm[key]
            val = data.obsm[key]
            if scipy.sparse.issparse(val):
                val = val.toarray()
                test_val = test_val.toarray()
            np.testing.assert_array_equal(test_val, val)


    # '--run-approximated-louvain',
    # '--run-louvain',
    # '--run-fle',
    # '--run-net-fle'
    # '--plot-hvg',


    def test_de_analysis(self):
        subprocess.run(
            ['scCloud', 'de_analysis', os.path.join('output', 'test_3k_pbmc.h5ad'), 'out_marker.de.xlsx', '--fisher',
             '--mwu', '--roc', '--labels', 'leiden_labels'])
        test_output = pd.read_excel('out_marker.de.xlsx', sheet_name=None)
        output = pd.read_excel(os.path.join('output', 'out_marker.de.xlsx'), sheet_name=None)
        for key in output:
            pd.testing.assert_frame_equal(test_output[key], output[key])
