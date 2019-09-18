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

	read_input
	write_output
	aggregate_matrices

Preprocess
~~~~~~~~~~


.. autosummary::
	:toctree: .

	qc_metrics
	get_filter_stats
	filter_data
	select_features
	log_norm
	highly_variable_features
	pca


Batch Correction
~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	set_group_attribute
	correct_batch

Nearest Neighbors
~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	neighbors
	calc_kBET
	calc_kSIM

Diffusion Map
~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	diffmap
	reduce_diffmap_to_3d
	calc_pseudotime
	infer_path


Cluster algorithms
~~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	louvain
	leiden
	spectral_louvain
	spectral_leiden

Visualization Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	tsne
	fitsne
	umap
	fle
	net_tsne
	net_fitsne
	net_umap
	net_fle

Differential Expression Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	de_analysis
	markers
	find_markers
	write_results_to_excel

Write single-cell-portal-formatted outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. module:: sccloud.tools
.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	tools.run_scp_output

Annotate clusters:
------------------

.. module:: sccloud
.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	infer_cell_types
	annotate

Plotting
--------

Interactive Plots
~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	embedding
	composition_plot
	variable_feature_plot
	heatmap
	dotplot

Quality Control Plots
~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	violin
	scatter
	scatter_matrix

Demultiplexing
--------------

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	estimate_background_probs
	demultiplex


Miscellaneous
-------------

.. currentmodule:: sccloud

.. autosummary::
	:toctree: .

	search_genes
	search_de_genes
