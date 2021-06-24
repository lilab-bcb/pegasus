1.4.0 :small:`June 24, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Add ``nmf`` and ``integrative_nmf`` functions to compute NMF and iNMF using ``nmf-torch``; ``integrative_nmf`` supports quantile normalization proposed in the ``LIGER`` papers.

* Change ``qc_metrics`` and now all defaults are ``None``.

* Fix a bug in ``heatmap`` plot.

* Fix bugs in ``spectral_leiden`` and ``spectral_louvain``.

* Improve message reporting execution time spent each step.

* Improve human immune cell markers and auto cell type assignment for human immune cells

* Add ``find_outlier_clusters`` function to detect if any cluster is an outlier regarding one of the qc attributes (n_genes, n_counts, percent_mito) using MWU test.

* ``annotate`` function now adds annotation as a categorical variable and sort categories in natural order.

* Update mouse brain cell markers

* ``scatter`` function now plots all cells if ``attrs`` == ``None``; Add ``fix_corners`` option to fix the four corners when only a subset of cells is shown.

* Support ``UMAP`` v0.5+.

* Update doublet detection algorithm.
