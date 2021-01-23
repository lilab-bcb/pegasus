pegasus demuxEM -p 2 --generate-diagnostic-plots tests/data/hashing_citeseq/cb_cc_raw_gene_bc_matrices_h5.h5 tests/data/hashing_citeseq/cb_cell_hashing.csv tests/cb_cc
pegasus aggregate_matrix tests/data/sample_hashing_citeseq.csv tests/cb_cc_citeseq
pegasus cluster -p 2 --select-singlets --louvain --umap --citeseq --citeseq-umap --citeseq-umap-exclude Mouse_IgG1,Mouse_IgG2a,Mouse_IgG2b,Rat_IgG2b tests/cb_cc_citeseq.zarr.zip tests/citeseq_result
pegasus plot scatter --basis umap --attributes louvain_labels,assignment tests/citeseq_result.zarr.zip tests/citeseq_result.umap.pdf
pegasus plot scatter --basis citeseq_umap --attributes louvain_labels,assignment tests/citeseq_result.zarr.zip tests/citeseq_result.citeseq_umap.pdf
