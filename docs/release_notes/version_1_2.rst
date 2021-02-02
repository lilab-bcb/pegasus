1.2.1 :small:`TBD`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Make PCA more reproducible. No need to keep options for robust PCA calculation:

    * In ``pca`` function, remove argument ``robust``.
    * In ``infer_doublets`` function, remove argument ``robust``.
    * In **pegasus cluster** command, remove option ``--pca-robust``.

* Add control on number of parallel threads for OpenMP/BLAS.

* Remove function ``reduce_diffmap_to_3d``. In **pegasus cluster** command, remove option ``--diffmap-to-3d``.

* Enhance plotting functions.

* Bug fix.

1.2.0 :small:`December 25, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* tSNE support: 
        
    * ``tsne`` function in API: Use `FIt-SNE <https://github.com/KlugerLab/FIt-SNE>`_ for tSNE embedding calculation. No longer support MulticoreTSNE.
    * Determine ``learning_rate`` argument in ``tsne`` more dynamically. ([Belkina19]_, [Kobak19]_)
    * By default, use PCA embedding for initialization in ``tsne``. ([Kobak19]_)
    * Remove ``net_tsne`` and ``fitsne`` functions from API.
    * Remove ``--net-tsne`` and ``--fitsne`` options from **pegasus cluster** command.

* Add multimodal support on RNA and CITE-Seq data back: ``--citeseq``, ``--citeseq-umap``, and ``--citeseq-umap-exclude`` in **pegasus cluster** command.

* Doublet detection:

    * Add automated doublet cutoff inference to ``infer_doublets`` function in API. ([Li20-2]_)
    * Expose doublet detection to command-line tool: ``--infer-doublets``, ``--expected-doublet-rate``, and ``--dbl-cluster-attr`` in **pegasus cluster** command.
    * Add doublet detection tutorial.

* Allow multiple marker files used in cell type annotation: ``annotate`` function in API; ``--markers`` option in **pegasus annotate_cluster** command.

* Rename *pc_regress_out* function in API to ``regress_out``.

* Update the regress out tutorial.

* Bug fix.
