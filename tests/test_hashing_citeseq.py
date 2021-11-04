import os
import glob
import unittest

import numpy as np
import pandas as pd
import pegasus as pg


class TestPipeline(unittest.TestCase):
    def test_demux(self):
        data = pg.read_input("tests/cb_cc_demux.zarr.zip")
        self.assertEqual(data.shape, (737280, 33694), "Demux data shape differs!")
        self.assertIn('demux_type', data.obs.columns, "Demux type is lost!")
        self.assertIn('assignment', data.obs.columns, "Cell assignment is lost!")
        f_list = glob.glob("tests/cb_cc.*.pdf")
        self.assertEqual(len(f_list), 4, "Demux diagnosis plots are missing!")
        self.assertIn('cb_cc.out.demuxEM.zarr.zip', os.listdir('tests'), "Demultiplexed RNA matrix is lost!")

    def test_citeseq(self):
        data = pg.read_input("tests/cb_cc_citeseq.zarr.zip")
        self.assertSetEqual(set(data.list_data()), set(['GRCh38-citeseq', 'GRCh38-rna']), "Some modality is missing!")
        self.assertIn('demux_type', data.obs.columns, "Demux type is lost!")
        self.assertIn('assignment', data.obs.columns, "Cell assignment is lost!")
        self.assertEqual(data.shape, (737280, 33694), "RNA data shape differs!")
        data.select_data('GRCh38-citeseq')
        self.assertEqual(data.shape, (578353, 31), "CITE-Seq data shape differs!")

    def test_clustering(self):
        data = pg.read_input("tests/citeseq_result.zarr.zip")
        self.assertSetEqual(set(data.list_data()), set(['GRCh38-citeseq', 'GRCh38-rna']), "Some modality is missing!")
        n_rna_cells = data.shape[0]
        self.assertNotIn('demux_type', data.obs.columns, "Demux type is not removed!")
        self.assertEqual(data.obs['assignment'].cat.categories.size, 7, "Not all cells are demultiplexed singlets!")
        self.assertIn('X_citeseq', data.obsm.keys(), "CITE-Seq coordinates are lost!")
        self.assertEqual(data.obsm['X_citeseq_umap'].shape[1], data.obsm['X_umap'].shape[1], "Some of UMAP embeddings is lost!")
        data.select_data('GRCh38-citeseq')
        n_citeseq_cells = data.shape[0]
        self.assertEqual(n_rna_cells, n_citeseq_cells, "Two modalities have inconsistent number of cells!")

    def test_plot(self):
        self.assertIn('citeseq_result.citeseq_umap.pdf', os.listdir('tests'), "CITE-Seq UMAP plot is lost!")
        self.assertIn('citeseq_result.umap.pdf', os.listdir('tests'), "RNA UMAP plot is lost!")

if __name__ == "__main__":
    unittest.main()
