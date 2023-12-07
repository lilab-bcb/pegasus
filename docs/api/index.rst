.. automodule:: pegasus
    :noindex:

API
===

*Pegasus* can also be used as a python package. Import pegasus by::

    import pegasus as pg

Read and Write
---------------

.. autosummary::
    :toctree: .

    read_input
    write_output
    aggregate_matrices

Analysis Tools
---------------

Preprocess
~~~~~~~~~~

.. autosummary::
    :toctree: .

    qc_metrics
    get_filter_stats
    filter_data
    identify_robust_genes
    log_norm
    log1p
    normalize
    arcsinh
    highly_variable_features
    select_features
    pca
    nmf
    regress_out
    calculate_z_score


Batch Correction
~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    run_harmony
    run_scanorama
    integrative_nmf
    run_scvi

Nearest Neighbors
~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    neighbors
    get_neighbors
    calc_kBET
    calc_kSIM

Diffusion Map
~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    diffmap
    calc_pseudotime
    infer_path


Cluster Algorithms
~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    cluster
    louvain
    leiden
    split_one_cluster
    spectral_louvain
    spectral_leiden

Visualization Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    tsne
    umap
    fle
    net_umap
    net_fle

Doublet Detection
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    infer_doublets
    mark_doublets

Gene Module Score
~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    calc_signature_score

Label Transfer
~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: .

    train_scarches_scanvi
    predict_scarches_scanvi

Differential Expression and Gene Set Enrichment Analysis
--------------------------------------------------------

.. autosummary::
    :toctree: .

    de_analysis
    markers
    write_results_to_excel
    fgsea

Annotate clusters
-------------------

.. autosummary::
    :toctree: .

    infer_cell_types
    annotate

Plotting
--------

.. autosummary::
    :toctree: .

    scatter
    scatter_groups
    spatial
    compo_plot
    violin
    heatmap
    dotplot
    dendrogram
    hvfplot
    qcviolin
    volcano
    rank_plot
    ridgeplot
    wordcloud
    plot_gsea
    elbowplot

Pseudo-bulk analysis
---------------------

.. autosummary::
    :toctree: .

    pseudobulk
    deseq2
    pseudo.markers
    pseudo.write_results_to_excel
    pseudo.volcano

Demultiplexing
---------------

.. autosummary::
    :toctree: .

    estimate_background_probs
    demultiplex
    attach_demux_results

Miscellaneous
-------------

.. autosummary::
    :toctree: .

    search_genes
    search_de_genes
    find_outlier_clusters
    find_markers
