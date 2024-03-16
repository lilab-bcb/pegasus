1.9.1 :small:`March 16, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Feature**

* Add ``write_fgsea_results_to_excel`` function (see `documentation <api/pegasus.write_fgsea_results_to_excel.html>`_)

**Improvement**

* For ``dotplot`` function, add ``show_only_expressed`` parameter to decide whether the color intensity of dots are based on cells expressing the genes to show or all cells. By default, ``show_only_expressed=True``. (PR `292 <https://github.com/lilab-bcb/pegasus/pull/292>`_)
* Allow ``scatter`` and ``spatial`` functions to have a list of ``vmin`` and ``vmax`` values when plotting multiple features.
* For ``scatter``, ``spatial``, ``dotplot``, ``violin`` functions, when some features are not in the data, emit a warning message and continue with features existing in the data.
* In ``spatial`` function, add ``nrows`` and ``ncols`` to organize subplots.
* In ``calculate_z_score`` function, enforce the input count matrix to be dense or ``scipy.csr_matrix``.
* In ``plot_gsea`` function, add ``label_fontsize`` parametere to allow change font label size in GSEA plots.

1.9.0 :small:`January 19, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Feature and Improvement**

* ``calculate_z_score`` works with sparse count matrix. [PR `276 <https://github.com/lilab-bcb/pegasus/pull/276>`_ Thanks to `Jayaram Kancherla <https://github.com/jkanche>`_]
* Plotting functions (``scatter``, ``dotplot``, ``violin``, ``heatmap``) now give warnings on genes/attributes not existing in the data, and skip them in the plots.
* Improve ``heatmap``:

  * Add ``show_sample_name`` parameter for cases of pseudo-bulk data, nanoString DSP data, etc.
  * Use Scipy's linkage (``scipy.cluster.hierarchy.linkage``) for dendrograms to use its optimal ordering feature for better results (see ``groupby_optimal_ordering`` parameter).

* Update human lung and mouse immune markers used by ``infer_cell_types`` function.
* ``run_harmony`` can accept multiple attributes to be the batch key, by providing a list of attribute names to its ``batch`` parameter.
* Expose ``online_batch_size`` parameter in ``nmf`` and ``integrative_nmf`` functions.
