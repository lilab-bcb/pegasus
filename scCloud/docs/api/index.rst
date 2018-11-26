API
===

``scCloud`` can also be used as a python package. Import ``scCloud`` by::

	import scCloud

Tools:
------

**Aggregate channel-specific count matrices**

.. autosummary::
	:toctree: .

	tools.aggregate_10x_matrices

**Preprocess**

.. autosummary::
	:toctree: .

	tools.read_input
	tools.update_var_names
	tools.filter_data
	tools.log_norm
	tools.run_pca
	tools.run_rpca
	tools.get_anndata_for_subclustering
	tools.filter_cells_cite_seq

**Batch correction**

.. autosummary::
	:toctree: .

	tools.set_group_attribute
	tools.estimate_adjustment_matrices
	tools.filter_genes_dispersion
	tools.collect_variable_gene_matrix
	tools.correct_batch_effects

**Diffusion map**

.. autosummary::
	:toctree: .

	tools.run_diffmap
	tools.run_pseudotime_calculation

**Cluster algorithms**

.. autosummary::
	:toctree: .

	tools.run_louvain
	tools.run_kmeans
	tools.run_approximated_louvain

**Visualization algorithms**

.. autosummary::
	:toctree: .

	tools.run_tsne
	tools.run_fitsne
	tools.run_umap
	tools.run_force_directed_layout

**Differential expression analysis**

.. autosummary::
	:toctree: .

	tools.run_de_analysis

**Write single-cell-portal-formatted outputs**

.. autosummary::
	:toctree: .

	tools.run_scp_output

Annotate clusters:
------------------

.. autosummary::
	:toctree: .

	annotate_cluster.annotate_clusters
	annotate_cluster.annotate_anndata_object

Plotting:
---------

**Static plots**

.. autosummary::
	:toctree: .

	plotting.plot_composition
	plotting.plot_scatter
	plotting.plot_scatter_groups
	plotting.plot_scatter_genes
	plotting.plot_scatter_gene_groups
	plotting.plot_heatmap

**Interactive plots**

.. autosummary::
	:toctree: .

	plotting.scatter
	plotting.scatter_real
	plotting.scatter3d
	plotting.scatter3d_real

**Quality control plots**

.. autosummary::
	:toctree: .

	plotting.plot_qc_violin

DemuxEM
-------

.. autosummary::
	:toctree: .

	demuxEM.estimate_background_probs
	demuxEM.demultiplex
	demuxEM.down_sampling
	demuxEM.plot_adt_hist
	demuxEM.plot_rna_hist
	demuxEM.plot_bar
	demuxEM.plot_violin
	demuxEM.plot_heatmap
	demuxEM.plot_dataframe_bar
	demuxEM.plot_down_sampling 

CITE-Seq
--------

.. autosummary::
	:toctree: .

	cite_seq.merge_rna_and_adt_data

Miscellaneous:
--------------

.. autosummary::
	:toctree: .

	misc.search_genes
	misc.search_de_genes
