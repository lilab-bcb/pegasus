import os
import numpy as np
import pandas as pd
import scanpy.api as sc
import xlsxwriter
from collections import Counter

from . import batch_correction


transfer_gene_name = [(358, 'ENSG00000268991', 'FAM231C.2'), (921, 'ENSG00000278139', 'AL358075.4'), (2207, 'ENSG00000232995', 'RGS5.2'), (5847, 'ENSG00000282827', 'AC134772.2'), (5938, 'ENSG00000271858', 'CYB561D2.2'), (6087, 'ENSG00000241572', 'PRICKLE2-AS1.2'), (7213, 'ENSG00000249428', 'CFAP99.2'), (9596, 'ENSG00000280987', 'MATR3.2'), (9605, 'ENSG00000279686', 'AC142391.1'), (10277, 'ENSG00000282913', 'BLOC1S5.2'), (10867, 'ENSG00000124593', 'AL365205.1'), (11619, 'ENSG00000268592', 'RAET1E-AS1.2'), (13877, 'ENSG00000231963', 'AL662864.1'), (16117, 'ENSG00000225655', 'BX255923.1'), (16938, 'ENSG00000282955', 'RABL6.2'), (17241, 'ENSG00000265264', 'TIMM10B.2'), (18626, 'ENSG00000282682', 'C11orf71.2'), (18984, 'ENSG00000282883', 'AKR1C3.2'), (19226, 'ENSG00000150076', 'CCDC7.2'), (19346, 'ENSG00000264404', 'BX547991.1'), (21184, 'ENSG00000282031', 'TMBIM4.2'), (21230, 'ENSG00000257815', 'LINC01481.2'), (22033, 'ENSG00000228741', 'SPATA13.2'), (22037, 'ENSG00000281899', 'AL359736.3'), (22654, 'ENSG00000274827', 'LINC01297.2'), (23662, 'ENSG00000273259', 'AL049839.2'), (24019, 'ENSG00000211974', 'AC245369.1'), (26919, 'ENSG00000279257', 'C17orf100.2'), (26962, 'ENSG00000187838', 'PLSCR3'), (27137, 'ENSG00000255104', 'AC005324.4'), (27884, 'ENSG00000263715', 'LINC02210-CRHR1'), (28407, 'ENSG00000281844', 'FBF1.2'), (30440, 'ENSG00000283027', 'CAPS.2'), (32648, 'ENSG00000235271', 'LINC01422.2')]

def update_var_names(adata, genome):
	if adata.var_names[0].startswith(genome + "_"):
		n = len(genome) + 1
		adata.var['gene_ids'] = [x[n:] for x in adata.var['gene_ids']]
		adata.var_names = pd.Index([x[n:] for x in adata.var_names])

	gsyms = adata.var_names.values
	
	if genome == "GRCh38":
		for pos, gid, gsym in transfer_gene_name:
			assert adata.var.iloc[pos, 0] == gid
			gsyms[pos] = gsym
	else:	
		dup_ids = Counter()
		for i in range(gsyms.size):
			idn = dup_ids[gsyms[i]]
			dup_ids[gsyms[i]] += 1
			if idn > 0:
				gsyms[i] = gsyms[i] + ".{}".format(idn)
	
	adata.var_names = pd.Index(gsyms)



def cluster(input_name, output_name, **kwargs):
	output_name = os.path.abspath(output_name)
	output_dir = os.path.dirname(output_name)
	output_base = '.' + os.path.basename(output_name)

	sc.settings.verbosity = 3
	sc.settings.writedir = sc.settings.figdir = output_dir + '/'

	if not kwargs['input_h5ad']:
		adata = sc.read_10x_h5(input_name + '_10x.h5', kwargs['genome'])
	else:
		adata = sc.read(input_name + '.h5ad')
	update_var_names(adata, kwargs['genome'])

	adata.obs['Channel'] = ['-'.join(x.split('-')[:-2]) for x in adata.obs_names]

	adata.obs['n_genes'] = adata.X.getnnz(axis = 1)
	adata.obs['n_counts'] = adata.X.sum(axis = 1).A1
	mito_genes = [name for name in adata.var_names if name.startswith(kwargs['mito_prefix'])]
	adata.obs['percent_mito'] = adata[:, mito_genes].X.sum(axis=1).A1 / adata.obs['n_counts'].values

	if kwargs['diagnosis']:
		sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save=output_base)
		sc.pl.scatter(adata, x='n_counts', y='percent_mito', save=output_base + ".1")
		sc.pl.scatter(adata, x='n_counts', y='n_genes', save=output_base + ".2")
		print("Plotted QC figures.")

	# Filter cells	
	if kwargs['filtr_xlsx'] is not None:
		writer = pd.ExcelWriter(kwargs['filtr_xlsx'], engine='xlsxwriter')
		tot_c = adata.obs['Channel'].value_counts()
	
	obs_index = np.logical_and.reduce((adata.obs['n_genes'] >= kwargs['min_genes'], 
									   adata.obs['n_genes'] < kwargs['max_genes'],
									   adata.obs['percent_mito'] < kwargs['percent_mito']))
	adata._inplace_subset_obs(obs_index)

	if kwargs['filtr_xlsx'] is not None:
		kept_c = adata.obs['Channel'].value_counts().reindex(tot_c.index, fill_value = 0)
		df = pd.DataFrame({'Kept' : kept_c, 'Filt' : tot_c - kept_c, 'Total' : tot_c})
		df = df[['Kept', 'Filt', 'Total']]
		df.sort_values('Kept', inplace = True)
		df.to_excel(writer, sheet_name = "Cell filtration stats")

	# Filter genes
	adata.var['n_cells'] = adata.X.getnnz(axis = 0)
	adata.var['percent_cells'] = adata.var['n_cells'] / adata.shape[0]
	adata.var['robust'] = adata.var['percent_cells'] >= kwargs['percent_cells']

	if kwargs['filtr_xlsx'] is not None:
		idx = adata.var['robust'] == False
		df = pd.DataFrame({'n_cells': adata.var.loc[idx, 'n_cells'], 'percent_cells': adata.var.loc[idx, 'percent_cells']})
		df.sort_values('n_cells', ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Gene filtration stats")
		writer.save()
		print("Filtration results are written.")

	var_index = (adata.var['n_cells'] > 0).values
	adata._inplace_subset_var(var_index)
	print("After filteration, {nc} cells and {ng} genes are kept. Among {ng} genes, {nrb} genes are robust.".format(nc = adata.shape[0], ng = adata.shape[1], nrb = adata.var['robust'].sum()))

	# Obtain normalized and logged gene expression matrix
	batch_correction.normalization(adata, kwargs['norm_count'])
	adata.X = adata.X.log1p()

	# Augment key barcode annotations for batch correction
	if not kwargs['input_h5ad']:
		df = pd.read_csv(input_name + '.attr.csv', header = 0, index_col = 0, dtype = str)
		df = df.loc[adata.obs['Channel'], kwargs['attrs']]
		df.index = adata.obs.index
		adata.obs[kwargs['attrs']] = df

	if kwargs['correct_batch']:
		if len(kwargs['groupby']) == 0:
			adata.obs['Group'] = 'one_group'
		else:
			adata.obs['Group'] = adata.obs[kwargs['groupby']].apply(lambda x: '+'.join([str(y) for y in x]), axis = 1)

	# Calculate batch correction adjustment matrices
	if kwargs['correct_batch']:
		batch_correction.estimate_adjustment_matrices(adata)

	adata.write(output_name + '.h5ad')
	print("Written log normalized data.")

	# Select variable genes
	# ~ 10^3, for 10^5, will include majority of expressed genes
	# filter_result contains a boolean list (gene_subset) within robust genes
	filter_result = batch_correction.filter_genes_dispersion(adata, kwargs['correct_batch'])
	if kwargs['diagnosis']:
		sc.pl.filter_genes_dispersion(filter_result, save=output_base)
		print("Plotted variable gene plot.")

	adata_c = batch_correction.collect_variable_gene_matrix(adata, filter_result.gene_subset)


	keys = ['louvain_labels']
	if kwargs['key'] is not None:
		keys.append(kwargs['key'])

	# Correct batch effects
	if kwargs['correct_batch']:
		batch_correction.correct_batch_effects(adata_c)
	
	sc.pp.scale(adata_c, max_value=10) # if x > max_value, x = max_value; since we have many zeros, negative values will have much smaller magnitudes than positive values.

	# PCA analysis
	sc.tl.pca(adata_c)
	adata_c.obsm['X_pca'] *= -1
	
	if kwargs['diagnosis']:
		sc.pl.pca_scatter(adata_c, right_margin=0.2, save=output_base)
		sc.pl.pca_variance_ratio(adata_c, save=output_base)
		sc.pl.pca_loadings(adata_c, save=output_base)
		print("Plotted PCA plots.")

	adata_c.write(output_name + '_var.h5ad')
	print("Written pca results.")

	# tSNE analysis
	sc.tl.tsne(adata_c, n_jobs = kwargs['nthreads'])
	if kwargs['diagnosis']:
		sc.pl.tsne(adata_c, save=output_base)
		print("Plotted tSNE plot.")

	adata_c.write(output_name + '_var.h5ad')
	print("Written tsne results.")

	# Louvain clustering
	sc.tl.louvain(adata_c, n_neighbors = 50, resolution = kwargs['resolution'], n_jobs = kwargs['nthreads'])
	adata_c.obs['louvain_labels'] = [str(int(x) + 1) for x in adata_c.obs['louvain_groups']]

	adata_c.write(output_name + '_var.h5ad')
	if kwargs['output_loom']:
		adata_c.write_loom(output_name + '_var.loom')
	print("Written louvain cluster results.")

	# Copy clustering information to the big matrix
	adata.obs['louvain_groups'] = adata_c.obs['louvain_groups']
	adata.obs['louvain_labels'] = adata_c.obs['louvain_labels']
	adata.obsm['X_tsne'] = adata_c.obsm['X_tsne']
 
	adata.write(output_name + '.h5ad')
	if kwargs['output_loom']:
		adata.write_loom(output_name + '.loom')
	print("Written full data.")

	# Generate tSNE plots
	if kwargs['legend_on_data']:
		sc.pl.tsne(adata_c, color = keys, save = output_base, legend_loc = "on data")
	else:
		sc.pl.tsne(adata_c, color = keys, save = output_base, legend_fontsize = 10)
