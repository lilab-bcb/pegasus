1.4.5 :small:`January 24, 2022`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Make several dependencies optional to meet with different use cases.

1.4.4 :small:`October 22, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Use PegasusIO *v0.4.0+* for data manipulation.

* Add ``calculate_z_score`` function to calculate standardized z-scored count matrix.

1.4.3 :small:`July 25, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Allow ``run_harmony`` function to use GPU for computation.

1.4.2 :small:`July 19, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Bug fix for ``--output-h5ad`` and ``--citeseq`` options in ``pegasus cluster`` command.

1.4.1 :small:`July 17, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Add NMF-related options to ``pegasus cluster`` command.

* Add word cloud graph plotting feature to ``pegasus plot`` command.

* ``pegasus aggregate_matrix`` command now allow sample-specific filtration with parameters set in the input CSV-format sample sheet.

* Update doublet detection method: ``infer_doublets`` and ``mark_doublets`` functions.

* Bug fix.

1.4.0 :small:`June 24, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Add ``nmf`` and ``integrative_nmf`` functions to compute NMF and iNMF using `nmf-torch <https://pypi.org/project/nmf-torch/>`_ package; ``integrative_nmf`` supports quantile normalization proposed in the *LIGER* papers ([Welch19]_, [Gao21]_).

* Change the parameter defaults of function ``qc_metrics``: Now all defaults are ``None``, meaning not performing any filtration on cell barcodes.

* In **Annotate Clusters** API functions:

    * Improve human immune cell markers and auto cell type assignment for human immune cells. (``infer_cell_types`` function)

    * Update mouse brain cell markers (``infer_cell_types`` function)

    * ``annotate`` function now adds annotation as a categorical variable and sort categories in natural order.

* Add ``find_outlier_clusters`` function to detect if any cluster is an outlier regarding one of the qc attributes (n_genes, n_counts, percent_mito) using MWU test.

* In **Plotting** API functions:

    * ``scatter`` function now plots all cells if ``attrs`` == ``None``; Add ``fix_corners`` option to fix the four corners when only a subset of cells is shown.

    * Fix a bug in ``heatmap`` plotting function.

* Fix bugs in functions ``spectral_leiden`` and ``spectral_louvain``.

* Improvements:

    * Support *umap-learn* v0.5+. (``umap`` and ``net_umap`` functions)

    * Update doublet detection algorithm. (``infer_doublets`` function)

    * Improve message reporting execution time spent each step. (Pegasus command line tool)
