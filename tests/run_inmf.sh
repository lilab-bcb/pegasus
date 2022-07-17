pegasus aggregate_matrix tests/data/count_matrix.csv tests/aggr

if [ -f "tests/aggr.zarr.zip" ]; then
    pegasus cluster -p 2 --output-h5ad --output-loom --correct-batch-effect --correction-method inmf --louvain --umap tests/aggr.zarr.zip tests/inmf_result
fi
