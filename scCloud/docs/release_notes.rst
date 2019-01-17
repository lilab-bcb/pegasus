Version 0.9.0 `January 17, 2019`
--------------------------------

In 'scCloud cluster', changed '--output-seurat-compatible' to '--make-output-seurat-compatible'. scCloud will not generate output_name.seurat.h5ad. Instead, output_name.h5ad should be able to convert to a Seurat object directly. In the seurat object, raw.data slot refers to the filtered count data, data slot refers to the log-normalized expression data, and scale.data refers to the variable-gene-selected, scaled data.
In 'scCloud cluster', added '--min-umis' and '--max-umis' options to filter cells based on UMI counts.
In 'scCloud cluster', '--output-filtration-results' option does not require a spreadsheet name anymore. In addition, added more statistics such as median number of genes per cell in the spreadsheet.
In 'scCloud cluster', added '--plot-filtration-results' and '--plot-filtration-figsize' to support plotting filtration results.
Improved documentation on 'scCloud cluster' outputs.
Added 'scCloud parquet' command to transfer h5ad file into a parquet file for web-based interactive visualization.

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

