import os
import unittest

import scCloud.pipeline


class TestPipeline(unittest.TestCase):
    def tearDown(self):
        os.remove('test.h5')
        # os.remove('test.h5ad')

    def test_run(self):
        scCloud.tools.aggregate_10x_matrices('data/aggregate_test.csv', restrictions=None, attributes=['Version'],
                                             output_file='test.h5',
                                             google_cloud=False, select_singlets=False, ngene=None, is_dropseq=None)

        scCloud.pipeline.run_pipeline('test.h5', 'test',
                                      cite_seq=False,
                                      output_filtration_results=False,
                                      output_seurat_compatible=False,
                                      output_loom=False,
                                      batch_correction=True,
                                      filt_xlsx=None,
                                      processed=False,
                                      submat_to_dense=True,
                                      select_variable_genes=True,
                                      genome='mm10',
                                      pca_key='X_pca',
                                      n_jobs=1,
                                      diffmap_full_speed=True,
                                      subcluster=False,
                                      run_approx_louvain=True,
                                      group_attribute='Version',
                                      min_genes=500,
                                      max_genes=6000,
                                      mito_prefix="MT-",
                                      percent_mito=0.1,
                                      percent_cells=0.0005,
                                      norm_count=1.00E+05,
                                      random_state=0,
                                      nPC=50,
                                      nDC=50,
                                      diffmap_K=100,
                                      diffmap_alpha=0.5,
                                      run_louvain=True,
                                      louvain_resolution=1.3,
                                      run_kmeans=True,
                                      run_hdbscan=False,
                                      run_fle=False,
                                      pseudotime=None,
                                      kmeans_n_clusters=10,
                                      louvain_affinity='W_norm',
                                      approx_louvain_ninit=20,
                                      approx_louvain_nclusters=30,
                                      approx_louvain_resolution=1.3,
                                      run_tsne=True,
                                      tsne_perplexity=30,
                                      run_fitsne=False,
                                      run_umap=False,
                                      umap_on_diffmap=False)
