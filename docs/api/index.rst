.. module:: sccloud

.. automodule:: sccloud
    :noindex:

API
===

``sccloud`` can also be used as a python package. Import ``sccloud`` by::

	import sccloud as scc

Analysis Tools
--------------


Read and Write
~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	sccloud.read_input
	sccloud.write_output
	sccloud.aggregate_matrices

Preprocess
~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	sccloud.qc_metrics
	sccloud.get_filter_stats
	sccloud.filter_data
	sccloud.select_features
	sccloud.log_norm
	sccloud.highly_variable_features
	sccloud.pca


Batch Correction
~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	sccloud.set_group_attribute
	sccloud.correct_batch

Nearest Neighbors
~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	sccloud.neighbors
	sccloud.calc_kBET
	sccloud.calc_kSIM

Diffusion Map
~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	sccloud.diffmap
	sccloud.reduce_diffmap_to_3d
	sccloud.calc_pseudotime
	sccloud.infer_path


Cluster algorithms
~~~~~~~~~~~~~~~~~~


.. autosummary::
	:toctree: .

	sccloud.louvain
	sccloud.leiden
	sccloud.spectral_louvain
	sccloud.spectral_leiden

Visualization Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~


.. autosummary::
	:toctree: .

	sccloud.tsne
	sccloud.fitsne
	sccloud.umap
	sccloud.fle
	sccloud.net_tsne
	sccloud.net_fitsne
	sccloud.net_umap
	sccloud.net_fle

Differential Expression Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. autosummary::
	:toctree: .

	sccloud.de_analysis
	sccloud.markers
	sccloud.find_markers
	sccloud.write_results_to_excel

Write single-cell-portal-formatted outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. autosummary::
	:toctree: .

	sccloud.tools.run_scp_output

Annotate clusters:
------------------


.. autosummary::
	:toctree: .

	sccloud.infer_cell_types
	sccloud.annotate

Plotting
--------

Interactive Plots
~~~~~~~~~~~~~~~~~


.. autosummary::
	:toctree: .

	sccloud.embedding
	sccloud.composition_plot
	sccloud.variable_feature_plot
	sccloud.heatmap
	sccloud.dotplot

Quality Control Plots
~~~~~~~~~~~~~~~~~~~~~


.. autosummary::
	:toctree: .

	sccloud.violin
	sccloud.scatter
	sccloud.scatter_matrix

Demultiplexing
--------------


.. autosummary::
	:toctree: .

	sccloud.estimate_background_probs
	sccloud.demultiplex


Miscellaneous
-------------


.. autosummary::
	:toctree: .

	sccloud.search_genes
	sccloud.search_de_genes
