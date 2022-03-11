1.5.0 :small:`March 9, 2022`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Features**

* Spatial data analysis:

    * Enable ``pegasus.read_input`` function to load 10x Visium data: set ``file_type="visium"`` option.
    * Add ``pegasus.spatial`` function to generate spatial plot for 10x Visium data.
    * Add *Spatial Analysis Tutorial* in `Tutorials <tutorials.html>`_.

* Pseudobulk analysis: see `summary <api/index.html#pseudo-bulk-analysis>`_

    * Add ``pegasus.pseudobulk`` function to generate pseudobulk matrix.
    * Add ``pegasus.deseq2`` function to perform pseudobulk differential expression (DE) analysis, which is a Python wrapper of `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_.

        * Requires `rpy2 <https://rpy2.github.io/doc.html>`_ and the original *DESeq2* R package installed.
    * Add ``pegasus.pseudo.markers``, ``pegasus.pseudo.write_results_to_excel`` and ``pegasus.pseudo.volcano`` functions for processing pseudobulk DE results.
    * Add *Pseudobulk Analysis Tutorial* in `Tutorials <tutorials.html>`_.

* Add ``pegasus.fgsea`` function to perform Gene Set Enrichment Analysis (GSEA) on DE results and plotting, which is a Python wrapper of `fgsea <http://bioconductor.org/packages/release/bioc/html/fgsea.html>`_.

    * Requires `rpy2 <https://rpy2.github.io/doc.html>`_ and the original *fgsea* R package installed.

**API Changes**

* Function ``correct_batch``, which implements the L/S adjustment batch correction method, is obsolete.
  We recommend using ``run_harmony`` instead, which is also the default of ``--correct-batch-effect`` option in ``pegasus cluster`` command.

* ``pegasus.highly_variable_features`` allows specify custom attribute key for batches (``batch`` option), and thus remove ``consider_batch`` option.
  To select HVGs without considering batch effects, simply use the default, or equivalently use ``batch=None`` option.

* Add ``dist`` option to ``pegasus.neighbors`` function to allow use distance other than L2. (Contribution by `hoondy <https://github.com/hoondy>`_ in `PR 233 <https://github.com/lilab-bcb/pegasus/pull/233>`_)

    * Available options: ``l2`` for L2 (default), ``ip`` for inner product, and ``cosine`` for cosine similarity.

* The kNN graph returned by ``pegasus.neighbors`` function is now stored in ``obsm`` field of the data object, no longer in ``uns`` field.
  Moreover, the kNN affinity matrix is stored in ``obsp`` field.

**Improvements**

* Adjust ``pegasus.write_output`` function to work with `Zarr <https://zarr.readthedocs.io>`_ v2.11.0+.
