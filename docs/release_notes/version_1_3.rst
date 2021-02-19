1.3.0 :small:`February 2, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Make PCA more reproducible. No need to keep options for robust PCA calculation:

    * In ``pca`` function, remove argument ``robust``.
    * In ``infer_doublets`` function, remove argument ``robust``.
    * In **pegasus cluster** command, remove option ``--pca-robust``.

* Add control on number of parallel threads for OpenMP/BLAS.

    * Now n_jobs = -1 refers to use all physical CPU cores instead of logcal CPU cores.

* Remove function ``reduce_diffmap_to_3d``. In **pegasus cluster** command, remove option ``--diffmap-to-3d``.

* Enhance compo_plot and dotplot functions' usability.

* Bug fix.
