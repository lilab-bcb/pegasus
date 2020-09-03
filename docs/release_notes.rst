Version 1.0.0 `September 1, 2020`
--------------------------------

* New features:

    * Use ``zarr`` file format to handle data, which has a better I/O performance in general.
    * Multi-modality support:

        * Data are manipulated in Multi-modal structure in memory.
        * Support focus analysis on Unimodal data, and appending other Unimodal data to it. (``--focus`` and ``--append`` options in ``cluster`` command)
    * Calculate signature / gene module scores. (`calc_signature_score <api/pegasus.calc_signature_score.html>`_)
    * `Doublet detection <api/index.html#doublet-detection>`_ based on `Scrublet <https://github.com/AllonKleinLab/scrublet>`_: ``run_scrublet``, ``infer_doublets``, and ``mark_singlets``.
    * Principal-Component-level regress out. (`pc_regress_out <api/pegasus.pc_regress_out.html>`_)
    * Batch correction using `Scanorama <https://github.com/brianhie/scanorama>`_. (`run_scanorama <api/pegasus.run_scanorama.html>`_)
    * Allow DE analysis with sample attribute as condition. (Set ``condition`` argument in `de_analysis <api/pegasus.de_analysis.html>`_)
    * Use static plots to show results (see `Plotting <api/index.html#plotting>`_):

        * Provide static plots: composition plot, embedding plot (e.g. tSNE, UMAP, FLE, etc.), dot plot, feature plot, and volcano plot;
        * Add more gene-specific plots: dendrogram, heatmap, violin plot, quality-control violin, HVF plot.
* Deprecations:

    * No longer support ``h5sc`` file format, which was the output format of ``aggregate_matrix`` command in Pegasus version ``0.x``.
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

Version 0.17.2 `June 26, 2020`
--------------------------------

* Make Pegasus compatible with *umap-learn* v0.4+.
* Use *louvain* 0.7+ for Louvain clustering.
* Update tutorial.

Version 0.17.1 `April 6, 2020`
--------------------------------

* Improve pegasus command-line tool log.
* Add human lung markers.
* Improve log-normalization speed.
* Provide robust version of PCA calculation as an option.
* Add signature score calculation API.
* Fix bugs.

Version 0.17.0 `March 10, 2020`
--------------------------------

* Support *anndata* 0.7 and *pandas* 1.0.

* Better ``loom`` format output writing function.

* Bug fix on ``mtx`` format output writing function.

* Update human immune cell markers.

* Improve ``pegasus scp_output`` command.

Version 0.16.11 `February 28, 2020`
------------------------------------

* Add ``--remap-singlets`` and ``--subset-singlets`` options to 'cluster' command.

* Allow reading ``loom`` file with user-specified batch key and black list.

Version 0.16.9 `February 17, 2020`
-----------------------------------

Allow reading ``h5ad`` file with user-specified batch key.

Version 0.16.8 `January 30, 2020`
-----------------------------------

Allow input annotated ``loom`` file.

Version 0.16.7 `January 28, 2020`
-----------------------------------

Allow input ``mtx`` files of more filename formats.

Version 0.16.5 `January 23, 2020`
-----------------------------------

Add Harmony algorithm for data integration.

Version 0.16.3 `December 17, 2019`
-----------------------------------

Add support for loading mtx files generated from BUStools.

Version 0.16.2 `December 8, 2019`
-----------------------------------

Fix bug in 'subcluster' command.

Version 0.16.1 `December 4, 2019`
-----------------------------------

Fix one bug in clustering pipeline.

Version 0.16.0 `December 3, 2019`
-----------------------------------

* Change options in 'aggregate_matrix' command: remove '--google-cloud', add '--default-reference'.

* Fix bug in '--annotation' option of 'annotate_cluster' command.

* Fix bug in 'net_fle' function with 3-dimension coordinates.

* Use **fisher** package version 0.1.9 or above, as modifications in our forked **fisher-modified** package has been merged into it.

Version 0.15.0 `October 2, 2019`
-----------------------------------

Rename package to *PegasusPy*, with module name *pegasus*.

Version 0.14.0 `September 17, 2019`
-----------------------------------

Provide Python API for interactive analysis.

Version 0.10.0 `January 31, 2019`
---------------------------------

Added 'find_markers' command to find markers using LightGBM.

Improved file loading speed and enabled the parsing of channels from barcode strings for cellranger aggregated h5 files.

Version 0.9.0 `January 17, 2019`
--------------------------------

In 'cluster' command, changed '--output-seurat-compatible' to '--make-output-seurat-compatible'. Do not generate output_name.seurat.h5ad.
Instead, output_name.h5ad should be able to convert to a Seurat object directly. In the seurat object, raw.data slot refers to the filtered
count data, data slot refers to the log-normalized expression data, and scale.data refers to the variable-gene-selected, scaled data.

In 'cluster' command, added '--min-umis' and '--max-umis' options to filter cells based on UMI counts.

In 'cluster' command, '--output-filtration-results' option does not require a spreadsheet name anymore. In addition, added more statistics such as median number of genes per cell in the spreadsheet.

In 'cluster' command, added '--plot-filtration-results' and '--plot-filtration-figsize' to support plotting filtration results.
Improved documentation on 'cluster command' outputs.

Added 'parquet' command to transfer h5ad file into a parquet file for web-based interactive visualization.

Version 0.8.0 `November 26, 2018`
---------------------------------

Added support for checking index collision for CITE-Seq/hashing experiments.

Version 0.7.0 `October 26, 2018`
--------------------------------

Added support for CITE-Seq analysis.

Version 0.6.0 `October 23, 2018`
--------------------------------

Renamed scrtools to scCloud.

Added demuxEM module for cell/nuclei-hashing.

Version 0.5.0 `August 21, 2018`
-------------------------------

Fixed a problem related AnnData.

Added support for BigQuery.

Version 0.4.0 `August 2, 2018`
------------------------------

Added mouse brain markers.

Allow aggregate matrix to take 'Sample' as attribute.

Version 0.3.0 `June 26, 2018`
-----------------------------

scrtools supports fast preprocessing, batch-correction, dimension reduction, graph-based clustering, diffusion maps, force-directed layouts, and differential expression analysis, annotate clusters, and plottings.
