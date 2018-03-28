#!/usr/bin/env python

from sys import argv, exit
import argparse
import numpy as np
import pandas as pd
import scanpy.api as sc
import batch_correction
import xlsxwriter
import scanpy_utils


parser = argparse.ArgumentParser(description = "Run scanpy to obtain clusters.", formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input", help = "Input hdf5 file in 10x format.", metavar = "input")
parser.add_argument("output_src", help = "Output unique string.", metavar = "output_src")

parser.add_argument("--output_dir", help = "Output directory.", metavar = "<directory>", default = "/ahg/regevdata/projects/Manton/exp")
parser.add_argument("-p", "--threads", help = "Number of threads.", metavar = "<value>", type = int, default = 1)
parser.add_argument("--genome", help = "Genome.", metavar = "<genome>", default = "GRCh38")

parser.add_argument("--correct-batch-effect", help = "Correct batch effect by channel ID.", action = "store_true")
parser.add_argument("--batch-id-func", help = "Function to extract batch IDs.", metavar = "<lambda function>", default = None)
parser.add_argument("--factor-id-func", help = "Function to extract factor IDs.", metavar = "<lambda function>", default = None)
parser.add_argument("--subcluster-analysis", help = "Perform subcluster analysis.", metavar = "<annotation> <comma_separated_cluster_list>", nargs = 2)
parser.add_argument("--customize", nargs = 2, help = "hue by customized function using sample name", metavar = "<name> <function>")

parser.add_argument("--output-filtration-results", help = "Output number of cells and genes filtered to a spreadsheet.", metavar = "<file>", default = None)

parser.add_argument("--legend-on-data", help = "Put legends on data.", action = "store_true")

parser.add_argument("--min-genes", help = "Minimum number of genes to keep a cell.", metavar = "<int>", type = int, default = 500)
parser.add_argument("--max-genes", help = "Maximum number of genes (boundary excluded) to keep a cell.", metavar = "<int>", type = int, default = 6000)
parser.add_argument("--percent-mito", help = "Maximum percentage of mitochondrial genes (boundary excluded) to keep a cell.", metavar = "<float>", type = float, default = 0.1)
parser.add_argument("--min-cells", help = "Minimum number of cells to keep a gene.", metavar = "<int>", type = int, default = 20)
parser.add_argument("--min-umis", help = "Minimum number of UMIs to keep a gene.", metavar = "<int>", type = int, default = 30)

parser.add_argument("--resolution", help = "Resolution parameter for louvain clustering.", metavar = "<float>", type = float, default = 1.3)

parser.add_argument("--no-per-cell-normalization", dest = "norm", help = "Do not normalize per cell before taking the log.", action = "store_false")
parser.add_argument("--counts-per-cell-after", help = "Total count per cell after normalization.", metavar = "<float>", type = float, default = 1e5)

args = parser.parse_args()



sc.settings.verbosity = 3
sc.settings.writedir = sc.settings.figdir = args.output_dir + "/"

out_str = ".{0}".format(args.output_src)
num_threads = int(args.threads)

if args.subcluster_analysis is not None:
	adata = sc.read(args.input)
	idx = np.isin(adata.obs[args.subcluster_analysis[0]], args.subcluster_analysis[1].split(','))
	adata = adata[idx]
	adata = sc.AnnData(adata.X, adata.obs['percent_mito'], adata.var['gene_ids'])
	print("{0} cells are selected.".format(adata.shape[0]))

	if args.correct_batch_effect:
		batch_correction.update_batch_id(adata, args.batch_id_func)
		batch_correction.update_factor_id(adata, args.factor_id_func)
	# ~ 10^3, for 10^5, will include majority of expressed genes
	filter_result = batch_correction.filter_genes_dispersion(adata, True, args.correct_batch_effect)
	sc.pl.filter_genes_dispersion(filter_result, save=out_str)
else:
	adata = sc.read_10x_h5(args.input, args.genome)
	scanpy_utils.update_var_names(adata, args.genome)

	mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
	adata.obs['n_genes'] = np.sum(adata.X > 0, axis = 1).A1
	adata.var['n_counts'] = np.sum(adata.X, axis = 0).A1

	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save=out_str)
	sc.pl.scatter(adata, x='n_counts', y='percent_mito', save=out_str + ".1")
	sc.pl.scatter(adata, x='n_counts', y='n_genes', save=out_str + ".2")
	print("Plotted preprocessing figures.")

	# update batch IDs
	if args.correct_batch_effect or args.output_filtration_results is not None:
		batch_correction.update_batch_id(adata, args.batch_id_func)
		batch_correction.update_factor_id(adata, args.factor_id_func)

	# Filter cells	
	if args.output_filtration_results is not None:
		writer = pd.ExcelWriter(args.output_filtration_results, engine='xlsxwriter')
		tot_c = adata.obs['batch'].value_counts()
	
	sc.pp.filter_cells(adata, min_genes = args.min_genes)
	adata = adata[adata.obs['n_genes'] < args.max_genes, :]
	adata = adata[adata.obs['percent_mito'] < args.percent_mito, :]

	if args.output_filtration_results is not None:
		kept_c = adata.obs['batch'].value_counts().reindex(tot_c.index, fill_value = 0)
		df = pd.DataFrame({'Kept' : kept_c, 'Filt' : tot_c - kept_c, 'Total' : tot_c})
		df = df[['Kept', 'Filt', 'Total']]
		df.sort_values('Kept', inplace = True)
		df.to_excel(writer, sheet_name = "Cell filtration stats")

	# Filter genes
	if args.output_filtration_results is not None:
		df = pd.DataFrame({'n_cells': (adata.X > 0).sum(axis = 0).A1, 'n_counts' : adata.var['n_counts']})

	sc.pp.filter_genes(adata, min_cells = args.min_cells)
	adata = adata[:, adata.var['n_counts'] >= args.min_umis]

	if args.output_filtration_results is not None:
		df = df.loc[df.index.difference(adata.var_names)] # get only filtered genes
		df.sort_values('n_cells', ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Gene filtration stats")
		writer.save()
		print("Filtration results are written.")

	print("After filteration, {0} cells are kept.".format(adata.shape[0]))

	# Normalize and select variable genes
	if args.norm:
		sc.pp.normalize_per_cell(adata, counts_per_cell_after = args.counts_per_cell_after) # normalize such that each cell has total count equal to median count across samples

	# ~ 10^3, for 10^5, will include majority of expressed genes
	filter_result = batch_correction.filter_genes_dispersion(adata, False, args.correct_batch_effect)
	sc.pl.filter_genes_dispersion(filter_result, save=out_str)
	sc.pp.log1p(adata)

keys = ["louvain_labels"]
if args.customize is not None:
	plot_key = args.customize[0]
	scanpy_utils.add_obs_column(adata, plot_key, args.customize[1])
	keys.append(plot_key)

sc.write(args.output_src, adata)
print("Written log normalized data.")

adata_c = adata.copy()
adata_c = adata_c[:, filter_result.gene_subset]

if args.correct_batch_effect:
	batch_correction.regress_out(adata_c)
sc.pp.scale(adata_c, max_value=10)

sc.tl.pca(adata_c)
adata_c.obsm['X_pca'] *= -1
sc.pl.pca_scatter(adata_c, right_margin=0.2, save=out_str)
sc.pl.pca_variance_ratio(adata_c, save=out_str)
sc.pl.pca_loadings(adata_c, save=out_str)

sc.write("{0}_var".format(args.output_src), adata_c)
print("Written pca results.")

sc.tl.tsne(adata_c, n_jobs = num_threads)
sc.pl.tsne(adata_c, save=out_str)
sc.write("{0}_var".format(args.output_src), adata_c)
print("Written tsne results.")

sc.tl.louvain(adata_c, n_neighbors=50, resolution = args.resolution, n_jobs = num_threads)
adata_c.obs['louvain_labels'] = [str(int(x) + 1) for x in adata_c.obs['louvain_groups']]

sc.write("{0}_var".format(args.output_src), adata_c)
print("Written louvain cluster results.")

adata.obs['louvain_groups'] = adata_c.obs['louvain_groups']
adata.obs['louvain_labels'] = adata_c.obs['louvain_labels']
adata.obsm['X_tsne'] = adata_c.obsm['X_tsne']
 
sc.write(args.output_src, adata)
print("Written full data.")

if args.legend_on_data:
	sc.pl.tsne(adata_c, color=keys, save=out_str, legend_loc="on data")
else:
	sc.pl.tsne(adata_c, color=keys, save=out_str, legend_fontsize=10)
