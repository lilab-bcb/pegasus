1.1.0 :small:`December 7, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Improve doublet detection in Scrublet-like way using automatic threshold selection strategy: ``infer_doublets``, and ``mark_doublets``. Remove *Scrublet* from dependency, and remove ``run_scrublet`` function.

* Enhance performance of log-normalization (``log_norm``) and signature score calculation (``calc_signature_score``).

* In ``pegasus cluster`` command, add ``--genome`` option to specify reference genome name for input data of ``dge``, ``csv``, or ``loom`` format.

* Update `Regress out tutorial <_static/tutorials/regress_out.html>`_.

* Add ``ridgeplot``.

* Improve plotting functions: ``heatmap``, and ``dendrogram``.

* Bug fix.
