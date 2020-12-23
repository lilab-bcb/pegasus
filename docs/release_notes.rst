Release Notes
---------------

.. role:: small

.. note::
    Also see the release notes of `PegasusIO <https://pegasusio.readthedocs.io/en/latest/release_notes.html>`_.

Version 1.2
~~~~~~~~~~~~~

.. include:: stable_release.rst

Version 1.1
~~~~~~~~~~~~~

1.1.0 :small:`December 7, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Improve doublet detection in Scrublet-like way using automatic threshold selection strategy: ``infer_doublets``, and ``mark_doublets``. Remove *Scrublet* from dependency, and remove ``run_scrublet`` function.

* Enhance performance of log-normalization (``log_norm``) and signature score calculation (``calc_signature_score``).

* In ``pegasus cluster`` command, add ``--genome`` option to specify reference genome name for input data of ``dge``, ``csv``, or ``loom`` format.

* Update `Regress out tutorial <_static/tutorials/regress_out.html>`_.

* Add ``ridgeplot``.

* Improve plotting functions: ``heatmap``, and ``dendrogram``.

* Bug fix.

Version 1.0
~~~~~~~~~~~~~~

1.0.0 :small:`September 22, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* New features:

    * Use ``zarr`` file format to handle data, which has a better I/O performance in general.
    * Multi-modality support:

        * Data are manipulated in Multi-modal structure in memory.
        * Support focus analysis on Unimodal data, and appending other Unimodal data to it. (``--focus`` and ``--append`` options in ``cluster`` command)
    * Calculate signature / gene module scores. (`calc_signature_score <api/pegasus.calc_signature_score.html>`_)
    * `Doublet detection <api/index.html#doublet-detection>`_ based on `Scrublet <https://github.com/AllonKleinLab/scrublet>`_: ``run_scrublet``, ``infer_doublets``, and ``mark_doublets``.
    * Principal-Component-level regress out. (`pc_regress_out <api/pegasus.pc_regress_out.html>`_)
    * Batch correction using `Scanorama <https://github.com/brianhie/scanorama>`_. (`run_scanorama <api/pegasus.run_scanorama.html>`_)
    * Allow DE analysis with sample attribute as condition. (Set ``condition`` argument in `de_analysis <api/pegasus.de_analysis.html>`_)
    * Use static plots to show results (see `Plotting <api/index.html#plotting>`_):

        * Provide static plots: composition plot, embedding plot (e.g. tSNE, UMAP, FLE, etc.), dot plot, feature plot, and volcano plot;
        * Add more gene-specific plots: dendrogram, heatmap, violin plot, quality-control violin, HVF plot.
* Deprecations:

    * No longer support ``h5sc`` file format, which was the output format of ``aggregate_matrix`` command in Pegasus version ``0.x``.
    * Remove ``net_fitsne`` function.
* API changes:

    * In cell quality-control, default percent of mitochondrial genes is changed from **10.0** to **20.0**. (``percent_mito`` argument in `qc_metrics <api/pegasus.qc_metrics.html>`_; ``--percent-mito`` option in ``cluster`` command)
    * Move gene quality-control out of ``filter_data`` function to be a separate step. (`identify_robust_genes <api/pegasus.identify_robust_genes.html>`_)
    * DE analysis now uses MWU test by default, not t test. (`de_analysis <api/pegasus.de_analysis.html>`_)
    * `infer_cell_types <api/pegasus.infer_cell_types.html>`_ uses MWU test as the default ``de_test``.
* Performance improvement:

    * Speed up MWU test in DE analysis, which is inspired by `Presto <https://github.com/immunogenomics/presto>`_.
    * Integrate Fisher's exact test via Cython in DE analysis to improve speed.
* Other highlights:

    * Make I/O and count matrices aggregation a dedicated package `PegasusIO <https://pegasusio.readthedocs.io>`_.
    * Tutorials:

        * Update `Analysis tutorial <_static/tutorials/pegasus_analysis.html>`_;
        * Add 3 more tutorials: one on `plotting library <_static/tutorials/plotting_tutorial.html>`_,
          one on `batch correction and data integration <_static/tutorials/batch_correction.html>`_,
          and one on `regress out <_static/tutorials/regress_out.html>`_.

Version 0.x
~~~~~~~~~~~~~~~

0.17.2 :small:`June 26, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Make Pegasus compatible with *umap-learn* v0.4+.
* Use *louvain* 0.7+ for Louvain clustering.
* Update tutorial.

0.17.1 :small:`April 6, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Improve pegasus command-line tool log.
* Add human lung markers.
* Improve log-normalization speed.
* Provide robust version of PCA calculation as an option.
* Add signature score calculation API.
* Fix bugs.

0.17.0 :small:`March 10, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Support *anndata* 0.7 and *pandas* 1.0.

* Better ``loom`` format output writing function.

* Bug fix on ``mtx`` format output writing function.

* Update human immune cell markers.

* Improve ``pegasus scp_output`` command.

0.16.11 :small:`February 28, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Add ``--remap-singlets`` and ``--subset-singlets`` options to 'cluster' command.

* Allow reading ``loom`` file with user-specified batch key and black list.

0.16.9 :small:`February 17, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Allow reading ``h5ad`` file with user-specified batch key.

0.16.8 :small:`January 30, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Allow input annotated ``loom`` file.

0.16.7 :small:`January 28, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Allow input ``mtx`` files of more filename formats.

0.16.5 :small:`January 23, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add Harmony algorithm for data integration.

0.16.3 :small:`December 17, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Add support for loading mtx files generated from BUStools.

0.16.2 :small:`December 8, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fix bug in 'subcluster' command.

0.16.1 :small:`December 4, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fix one bug in clustering pipeline.

0.16.0 :small:`December 3, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Change options in 'aggregate_matrix' command: remove '--google-cloud', add '--default-reference'.

* Fix bug in '--annotation' option of 'annotate_cluster' command.

* Fix bug in 'net_fle' function with 3-dimension coordinates.

* Use **fisher** package version 0.1.9 or above, as modifications in our forked **fisher-modified** package has been merged into it.

0.15.0 :small:`October 2, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rename package to *PegasusPy*, with module name *pegasus*.

0.14.0 :small:`September 17, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Provide Python API for interactive analysis.

0.10.0 :small:`January 31, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added 'find_markers' command to find markers using LightGBM.

Improved file loading speed and enabled the parsing of channels from barcode strings for cellranger aggregated h5 files.

0.9.0 :small:`January 17, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In 'cluster' command, changed '--output-seurat-compatible' to '--make-output-seurat-compatible'. Do not generate output_name.seurat.h5ad.
Instead, output_name.h5ad should be able to convert to a Seurat object directly. In the seurat object, raw.data slot refers to the filtered
count data, data slot refers to the log-normalized expression data, and scale.data refers to the variable-gene-selected, scaled data.

In 'cluster' command, added '--min-umis' and '--max-umis' options to filter cells based on UMI counts.

In 'cluster' command, '--output-filtration-results' option does not require a spreadsheet name anymore. In addition, added more statistics such as median number of genes per cell in the spreadsheet.

In 'cluster' command, added '--plot-filtration-results' and '--plot-filtration-figsize' to support plotting filtration results.
Improved documentation on 'cluster command' outputs.

Added 'parquet' command to transfer h5ad file into a parquet file for web-based interactive visualization.

0.8.0 :small:`November 26, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added support for checking index collision for CITE-Seq/hashing experiments.

0.7.0 :small:`October 26, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added support for CITE-Seq analysis.

0.6.0 :small:`October 23, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Renamed scrtools to scCloud.

Added demuxEM module for cell/nuclei-hashing.

0.5.0 :small:`August 21, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed a problem related AnnData.

Added support for BigQuery.

0.4.0 :small:`August 2, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added mouse brain markers.

Allow aggregate matrix to take 'Sample' as attribute.

0.3.0 :small:`June 26, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

scrtools supports fast preprocessing, batch-correction, dimension reduction, graph-based clustering, diffusion maps, force-directed layouts, and differential expression analysis, annotate clusters, and plottings.
