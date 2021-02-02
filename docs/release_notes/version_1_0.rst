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

    * Make I/O and count matrices aggregation a dedicated package `PegasusIO <https://pegasusio.readthedocs.io>`__.
    * Tutorials:

        * Update `Analysis tutorial <_static/tutorials/pegasus_analysis.html>`_;
        * Add 3 more tutorials: one on `plotting library <_static/tutorials/plotting_tutorial.html>`_,
          one on `batch correction and data integration <_static/tutorials/batch_correction.html>`_,
          and one on `regress out <_static/tutorials/regress_out.html>`_.
