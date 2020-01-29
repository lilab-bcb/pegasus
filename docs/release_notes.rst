Version 0.16.7 `January 28, 2020`
-----------------------------------

Allow input ``mtx`` files of more filename formats.

Version 0.16.5 `January 23, 2020`
-----------------------------------

Add Harmony algorithm for data integration.

Version 0.16.3 `December 17, 2019`
-----------------------------------

Add support for loading mtx files generated from BUStools.

Version 0.16.2 `December 8, 2019`
-----------------------------------

Fix bug in 'subcluster' command.

Version 0.16.1 `December 4, 2019`
-----------------------------------

Fix one bug in clustering pipeline.

Version 0.16.0 `December 3, 2019`
-----------------------------------

* Change options in 'aggregate_matrix' command: remove '--google-cloud', add '--default-reference'.

* Fix bug in '--annotation' option of 'annotate_cluster' command.

* Fix bug in 'net_fle' function with 3-dimension coordinates.

* Use **fisher** package version 0.1.9 or above, as modifications in our forked **fisher-modified** package has been merged into it.

Version 0.15.0 `October 2, 2019`
-----------------------------------

Rename package to *PegasusPy*, with module name *pegasus*.

Version 0.14.0 `September 17, 2019`
-----------------------------------

Provide Python API for interactive analysis.

Version 0.10.0 `January 31, 2019`
---------------------------------

Added 'find_markers' command to find markers using LightGBM.

Improved file loading speed and enabled the parsing of channels from barcode strings for cellranger aggregated h5 files.

Version 0.9.0 `January 17, 2019`
--------------------------------

In 'cluster' command, changed '--output-seurat-compatible' to '--make-output-seurat-compatible'. Do not generate output_name.seurat.h5ad.
Instead, output_name.h5ad should be able to convert to a Seurat object directly. In the seurat object, raw.data slot refers to the filtered
count data, data slot refers to the log-normalized expression data, and scale.data refers to the variable-gene-selected, scaled data.

In 'cluster' command, added '--min-umis' and '--max-umis' options to filter cells based on UMI counts.

In 'cluster' command, '--output-filtration-results' option does not require a spreadsheet name anymore. In addition, added more statistics such as median number of genes per cell in the spreadsheet.

In 'cluster' command, added '--plot-filtration-results' and '--plot-filtration-figsize' to support plotting filtration results.
Improved documentation on 'cluster command' outputs.

Added 'parquet' command to transfer h5ad file into a parquet file for web-based interactive visualization.

Version 0.8.0 `November 26, 2018`
---------------------------------

Added support for checking index collision for CITE-Seq/hashing experiments.

Version 0.7.0 `October 26, 2018`
--------------------------------

Added support for CITE-Seq analysis.

Version 0.6.0 `October 23, 2018`
--------------------------------

Renamed scrtools to scCloud.

Added demuxEM module for cell/nuclei-hashing.

Version 0.5.0 `August 21, 2018`
-------------------------------

Fixed a problem related AnnData.

Added support for BigQuery.

Version 0.4.0 `August 2, 2018`
------------------------------

Added mouse brain markers.

Allow aggregate matrix to take 'Sample' as attribute.

Version 0.3.0 `June 26, 2018`
-----------------------------

scrtools supports fast preprocessing, batch-correction, dimension reduction, graph-based clustering, diffusion maps, force-directed layouts, and differential expression analysis, annotate clusters, and plottings.
