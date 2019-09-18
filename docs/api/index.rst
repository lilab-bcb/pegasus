.. automodule:: sccloud
    :noindex:

.. module:: sccloud


API
===

``sccloud`` can also be used as a python package. Import ``sccloud`` by::

	import sccloud as scc

Analysis Tools
--------------

Read and Write
~~~~~~~~~~~~~~

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

.. autosummary::
	:toctree: .

	set_group_attribute
	correct_batch

Nearest Neighbors
~~~~~~~~~~~~~~~~~

.. autosummary::
	:toctree: .

	neighbors
	calc_kBET
	calc_kSIM

Diffusion Map
~~~~~~~~~~~~~

.. autosummary::
	:toctree: .

	diffmap
	reduce_diffmap_to_3d
	calc_pseudotime
	infer_path


Cluster algorithms
~~~~~~~~~~~~~~~~~~

.. autosummary::
	:toctree: .

	louvain
	leiden
	spectral_louvain
	spectral_leiden

Visualization Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~

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

.. autosummary::
	:toctree: .

	de_analysis
	markers
	find_markers
	write_results_to_excel

Annotate clusters:
------------------

.. autosummary::
	:toctree: .

	infer_cell_types
	annotate

Plotting
--------

Interactive Plots
~~~~~~~~~~~~~~~~~

.. autosummary::
	:toctree: .

	embedding
	composition_plot
	variable_feature_plot
	heatmap
	dotplot

Quality Control Plots
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
	:toctree: .

	violin
	scatter
	scatter_matrix

Demultiplexing
--------------

.. autosummary::
	:toctree: .

	estimate_background_probs
	demultiplex


Miscellaneous
-------------

.. autosummary::
	:toctree: .

	search_genes
	search_de_genes
