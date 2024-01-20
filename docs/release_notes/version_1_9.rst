1.9.0 :small:`January 19, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Feature and Improvement**

* ``calculate_z_score`` works with sparse count matrix. [PR `276 <https://github.com/lilab-bcb/pegasus/pull/276>`_ Thanks to `Jayaram Kancherla <https://github.com/jkanche>`_]
* Plotting functions (``scatter``, ``dotplot``, ``violin``, ``heatmap``) now give warnings on genes/attributes not existing in the data, and skip them in the plots.
* Improve ``heatmap``:

  * Add ``show_sample_name`` parameter for cases of pseudo-bulk data, nanoString DSP data, etc.
  * Use Scipy's linkage (``scipy.cluster.hierarchy.linkage``) for dendrograms to use its optimal ordering feature for better results (see ``groupby_optimal_ordering`` parameter).

* Update human lung and mouse immune markers used by ``infer_cell_types`` function.
* Expose ``online_batch_size`` parameter in ``nmf`` and ``integrative_nmf`` functions.
