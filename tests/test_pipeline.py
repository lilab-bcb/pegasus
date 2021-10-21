import os
import glob
import unittest

import numpy as np
import pandas as pd
import pegasus as pg


class TestPipeline(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestPipeline, self).__init__(*args, **kwargs)
        self.aggr_data = pg.read_input("tests/aggr.zarr.zip")
        self.data = pg.read_input("tests/result.zarr.zip")

    def test_aggregate(self):
        self.assertEqual(self.aggr_data.shape, (1723, 31053), "Shape differs!")
        self.assertEqual(self.data.obs['Channel'].cat.categories.size, 2, "Number of samples differs!")
        self.assertEqual(self.data._unidata.get_uid(), 'mm10-rna', "Genome and Modality info differs!")

    def test_qc(self):
        self.assertEqual(self.data.shape[0], 1043, "Number of cell barcodes differs!")
        f_list = glob.glob("tests/*.filt.*.pdf")
        self.assertEqual(len(f_list), 3, "Number of QC plots differs!")
        self.assertIn('df_qcplot', self.data.uns, "QC plot cache is lost!")

    def test_clustering(self):
        self.assertEqual(self.data.obsm['pca_harmony_knn_indices'].shape, (1043, 99), "KNN graph shape differs!")
        self.assertEqual(self.data.obsm['pca_harmony_knn_distances'].shape, (1043, 99), "KNN distance matrix shape differs!")
        self.assertIn('louvain_labels', self.data.obs.columns, "Louvain result is lost!")
        self.assertIn('leiden_labels', self.data.obs.columns, "Leiden result is lost!")

    def test_doublet_detection(self):
        self.assertIn('doublet_score', self.data.obs.columns, "Doublet score is lost!")
        self.assertEqual(self.data.obs['pred_dbl'].dtype, 'bool', "Type differs!")
        self.assertSetEqual(set(self.data.obs['demux_type'].cat.categories.tolist()), set(['doublet', 'singlet']), "Demux type names differ!")

    def test_signature_score(self):
        self.assertIn('G1/S', self.data.obs.columns, "G1/S signature score is lost!")
        self.assertIn('G2/M', self.data.obs.columns, "G2/M signature score is lost!")
        self.assertIn('cycling', self.data.obs.columns, "Max of G1/S and G2/M signature score is lost!")
        self.assertSetEqual(set(self.data.obs['predicted_phase'].cat.categories.tolist()), set(['G0', 'G1/S', 'G2/M']), "Predicted phases of cell cycling differ!")

    def test_embeddings(self):
        self.assertEqual(self.data.obsm['X_pca_harmony'].shape, (1043, 50), "Harmony PCA shape differs!")
        self.assertEqual(self.data.obsm['X_umap'].shape, (1043, 2), "UMAP shape differs!")
        self.assertEqual(self.data.obsm['X_tsne'].shape, (1043, 2), "tSNE shape differs!")
        self.assertEqual(self.data.obsm['X_fle'].shape, (1043, 2), "FLE shape differs!")

    def test_de_tests(self):
        df_de = pd.DataFrame(self.data.varm['de_res'], index=self.data.var_names)
        self.assertIn('1:auroc', df_de.columns)
        self.assertIn('1:mwu_qval', df_de.columns)
        self.assertIn('1:t_qval', df_de.columns)
        self.assertIn('1:fisher_qval', df_de.columns)
        self.assertIn('result.de.xlsx', os.listdir('tests'), "DE spread sheet is lost!")

    def test_annotation(self):
        self.assertIn('result.anno.txt', os.listdir('tests'), "Cell type annotation result is lost!")

    def test_plot(self):
        self.assertIn('result.compo.pdf', os.listdir('tests'), "Composition plot is lost!")
        self.assertIn('result.louvain_labels.umap.pdf', os.listdir('tests'), "UMAP plot is lost!")
        self.assertIn('result.leiden_labels.tsne.pdf', os.listdir('tests'), "tSNE plot is lost!")
        self.assertIn('result.louvain_labels.fle.pdf', os.listdir('tests'), 'FLE plot is lost!')

    def test_output(self):
        data_h5ad = pg.read_input("tests/result.mm10-rna.h5ad")
        self.assertEqual(self.data.shape, data_h5ad.shape, "H5AD format's shape is inconsistent!")

        data_loom = pg.read_input("tests/result.mm10-rna.loom")
        self.assertEqual(self.data.shape, data_loom.shape, "Loom format's shape is inconsistent!")

        self.assertIn('result.log', os.listdir('tests'), 'Clustering log is lost!')


if __name__ == "__main__":
    unittest.main()
