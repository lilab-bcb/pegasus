pegasus aggregate_matrix tests/data/count_matrix.csv tests/aggr
pegasus cluster -p 2 --min-genes 500 --max-genes 6000 --percent-mito 20.0 --output-filtration-results --output-h5ad --output-loom --plot-filtration-results --plot-hvf --correct-batch-effect --pca-robust --louvain --leiden --tsne --umap --fle --infer-doublets --dbl-cluster-attr louvain_labels tests/aggr.zarr.zip tests/result
pegasus de_analysis -p 2 --labels louvain_labels --t --fisher tests/result.zarr.zip tests/result.de.xlsx
pegasus annotate_cluster --markers mouse_immune,mouse_brain tests/result.zarr.zip tests/result.anno.txt
pegasus plot compo --groupby leiden_labels --condition Channel tests/result.zarr.zip tests/result.compo.pdf
pegasus plot scatter --basis umap --attributes louvain_labels,Channel tests/result.zarr.zip tests/result.louvain_labels.umap.pdf
pegasus plot scatter --basis tsne --attributes leiden_labels,Channel tests/result.zarr.zip tests/result.leiden_labels.tsne.pdf
pegasus plot scatter --basis fle --attributes louvain_labels,Channel tests/result.zarr.zip tests/result.louvain_labels.fle.pdf
