1.10.0 :small:`June 12, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Features**

* Add ``pegasus.pseudo.get_original_DE_result`` function to return the DE result as a Pandas DataFrame.
* Implement a new version of Dendrogram (PR `295  <https://github.com/lilab-bcb/pegasus/pull/295>`_):

  * It now uses Connection Specific Index (CSI) for the distance calculation.
  * A new function ``calc_dendrogram`` to calculate the linkage, and store it in ``data.uns`` field.
  * A new function ``plot_dendrogram`` to plot the dendrogram based on the linkage calculated.
* For ``deseq2`` function (PR `300 <https://github.com/lilab-bcb/pegasus/pull/300>`_, `304 <https://github.com/lilab-bcb/pegasus/pull/304>`_):

  * It now has 2 backends: ``pydeseq2`` by default, which uses Python package PyDESeq2; ``deseq2`` for R package DESeq2. Specify it in ``backend`` parameter.
  * Add ``compute_all`` parameter to decide if applying DE analysis to all count matrices of the input pseudobulk data object, or only for the default count matrix.
  * Add ``alpha`` parameter o allow user choose significance level for independent filtering to calculate adjusted p-values.
  * For ``pydeseq2`` backend, it accepts inference on multiple contrasts in the same model.
* Reorganize GSEA functions (PR `305 <https://github.com/lilab-bcb/pegasus/pull/305>`_):

  * Rename ``fgsea`` function to ``gsea``. The new ``gsea`` function accepts two methods: ``gseapy`` to use GSEAPy's prerank function, which is the default; ``fgsea`` to use R package fgsea.
  * Use parameter ``rank_key`` to specify which attribute in the DE result to be used as gene ranks/signatures.
  * The number of threads ``n_jobs`` is ``4`` by default.
  * Rename ``write_fgsea_results_to_excel`` function to ``write_gsea_results_to_excel``.

**Improvement**

* Add ``label_fontsize`` parameter to ``plot_gsea`` function to allow change label font size in GSEA plots.
* Improve ``scatter`` function:

  * Add ``aspect`` parameter. By default it's ``aspect=auto`` to enforce square plots; set to ``aspect=equal`` for the actual y to x axis ratio.
* Improve ``spatial`` function (PR `294 <https://github.com/lilab-bcb/pegasus/pull/294>`_):

  * Now accepts non-Visium data, which may not have spatial images stored in ``data.img`` field.
  * Add parameter ``aspect``. Default is ``aspect=equal`` to plot the actual y to x axis ratio; set to ``aspect=auto`` to enforce square plots.
  * Add parameter ``margin_percent`` to set the image margin to the 4 sides.
  * Add parameters ``restrictions``, ``show_background`` and ``palettes``, which have the same behaviors as in ``scatter`` function.
  * Add parameter ``y_flip`` with default ``True``. This is for the case where y coordinates start from top, which needs a flip. Set to ``False`` if the spatial y-coordinates start from bottom.
  * When plotting Visium data with no image background: set ``resolution=None`` and ``y_flip=True``.
* For ``pseudobulk`` function (` <(PR `300 <https://github.com/lilab-bcb/pegasus/pull/300>`_)>`_): rename parameter ``sample`` to ``groupby``, ``cluster`` to ``condition``; it returns a ``MultimodalData`` object, and does not write to the input data object.
* For ``pegasus.pseudo.volcano`` function, add ``rank_by`` parameter to rank genes by log2 fold change (``log2fc``) or -log10(p-value) (``neglog10p``). This would affect the top gene names shown in the plots.
* Bug fix in plotting functions to work with Matplotlib v3.9+.
* Add bimodality index to doublet detection (PR `307 <https://github.com/lilab-bcb/pegasus/pull/307>`_).
