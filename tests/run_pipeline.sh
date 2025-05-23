pegasus aggregate_matrix tests/data/count_matrix.csv tests/aggr

if [ -f "tests/aggr.zarr.zip" ]; then
    pegasus cluster -p 2 --min-genes 500 --max-genes 6000 --mito-prefix mt- --percent-mito 20.0 --output-filtration-results --output-h5ad --output-loom --plot-filtration-results --plot-hvf --exact-K --correct-batch-effect --nmf --leiden --umap --fle --calc-signature-scores cell_cycle_mouse tests/aggr.zarr.zip tests/result
fi

if [ -f "tests/result.zarr.zip" ]; then
    pegasus de_analysis -p 2 --labels leiden_labels --t --fisher tests/result.zarr.zip tests/result.de.xlsx
    pegasus annotate_cluster --markers mouse_immune,mouse_brain tests/result.zarr.zip tests/result.anno.txt
    pegasus plot compo --groupby leiden_labels --condition Channel tests/result.zarr.zip tests/result.compo.pdf
    pegasus plot scatter --basis umap --attributes leiden_labels,Channel tests/result.zarr.zip tests/result.leiden_labels.umap.pdf
    pegasus plot scatter --basis fle --attributes leiden_labels,Channel tests/result.zarr.zip tests/result.leiden_labels.fle.pdf
fi
