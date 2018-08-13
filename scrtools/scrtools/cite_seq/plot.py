import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages



def plot_barcode_hist(data, adt, out_file, dpi = 500):
	background = adt.X.sum(axis = 1).A1

	idx_df = data.obs_names.isin(adt.obs_names)
	if idx_df.sum() < data.shape[0]:
		nzero = data.shape[0] - idx_df.sum()
		print("Warning: {} cells do not have ADTs, percentage = {:.2f}%.".format(nzero, nzero * 100.0 / data.shape[0]))

	signal = adt[data.obs_names[idx_df],].X.sum(axis = 1).A1
	bins = np.logspace(0, np.log10(max(background)), 101)
	plt.hist(background, bins, alpha = 0.5, label = 'background', log = True)
	plt.hist(signal, bins, alpha = 0.5, label = 'signal', log = True)
	plt.legend(loc='upper right')
	ax = plt.gca()
	ax.set_xscale("log")
	ax.set_xlabel("Number of UMIs (log10 scale)")
	ax.set_ylabel("Number of barcodes (log10 scale)")
	plt.savefig(out_file, dpi = dpi)
	plt.close()



def plot_antibody_hist(adts, names, antibody, control, out_file, dpi = 500, figsize = None, format = None):
	fig, axes = plt.subplots(nrows = 2, ncols = len(adts), squeeze = False, figsize = figsize)
	plt.tight_layout(pad = 4)

	for i, adt in enumerate(adts):
		signal = adt[:, antibody].X.toarray()
		background = adt[:, control].X.toarray()
		bins = np.logspace(0, np.log10(max(signal.max(), background.max())), 101)

		ax = axes[0, i]
		ax.hist(background, bins, alpha = 0.5, label = control, log = True)
		ax.hist(signal, bins, alpha = 0.5, label = antibody, log = True)
		ax.legend(loc='upper right')
		ax.set_xscale("log")
		ax.set_xlabel("Number of UMIs (log10 scale)")
		ax.set_ylabel("Number of barcodes (log10 scale)")
		ax.set_title(names[i])

		idx = (signal > 0) | (background > 0)
		fc = (signal[idx] + 1.0) / (background[idx] + 1.0)
		bins = np.logspace(np.log10(fc.min()), np.log10(fc.max()), 101)

		ax = axes[1, i]
		ax.hist(fc, bins, label = "{} /\n{}".format(antibody, control), log = True)
		ax.legend(loc='upper right')
		ax.set_xscale("log")
		ax.set_xlabel("Fold change (log10 scale)")
		ax.set_ylabel("Number of barcodes (log10 scale)")
		ax.set_title("{:.2%}".format((fc >= 10.0).sum() / fc.size))

	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()

def plot_antibodies_hist(adts, names, df, out_file, figsize = None):
	with PdfPages(out_file) as pdf:
		for idx, row in df.iterrows():
			plot_antibody_hist(adts, names, row['antibody'], row['control'], pdf, figsize = figsize, format = "pdf")
			print("{} / {}".format(row['antibody'], row['control']))



def plot_gene_violin(data, gene_list, out_file, dpi = 500, figsize = (6, 4)):
	df = pd.DataFrame(data.X.toarray(), index = data.obs_names, columns = data.var_names)
	df['assignment'] = data.obs['assignment']

	gene_list = np.array(gene_list)
	gene_list = gene_list[np.isin(gene_list, df.columns)]

	nrows = int(np.sqrt(gene_list.size) - 1e-3) + 1
	ncols = int(gene_list.size / nrows - 1e-3) + 1

	fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize = figsize, squeeze = False, frameon = False)
	for i in range(nrows):
		for j in range(ncols):
			ax = axes[i][j]
			ax.grid(False)
			gid = i * ncols + j
			if gid >= gene_list.size:
				ax.set_frame_on(False)
				ax.set_xticks([])
				ax.set_yticks([])
				continue
			sns.violinplot(x = 'assignment', y = gene_list[gid], data = df, ax = ax)
			ax.set_ylabel('log TPM')
			ax.set_title(gene_list[gid])

	plt.tight_layout()
	plt.savefig(out_file, dpi = dpi)
	plt.close()



def plot_doublet_hist(adata, nuclei_type, out_file, attr = 'coarse_annotation', dpi = 500, format = None):
	fig, axes = plt.subplots(nrows = 2, ncols = 2, squeeze = False)
	plt.tight_layout(pad = 2)

	idx = np.isin(adata.obs[attr], nuclei_type)	
	bins = np.linspace(adata.obs.loc[idx, 'n_genes'].min(), adata.obs.loc[idx, 'n_genes'].max(), 101)
	
	ax = axes[0, 0]
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'z_unknown'), 'n_genes'], bins, alpha = 0.5, label = 'unknown')
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'z_doublet'), 'n_genes'], bins, alpha = 0.5, label = 'doublet')
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'hashtag'), 'n_genes'], bins, alpha = 0.5, label = 'singlet')
	ax.legend(loc='upper right')
	ax.set_xlabel("Number of genes")
	ax.set_ylabel("Number of nuclei")

	ax = axes[0, 1]
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'nohashtag'), 'n_genes'], bins, label = 'nohash')
	ax.legend(loc='upper right')
	ax.set_xlabel("Number of genes")
	ax.set_ylabel("Number of nuclei")

	bins = np.linspace(adata.obs.loc[idx, 'n_counts'].min(), adata.obs.loc[idx, 'n_counts'].max(), 101)
	
	ax = axes[1, 0]
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'z_unknown'), 'n_counts'], bins, alpha = 0.5, label = 'unknown')
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'z_doublet'), 'n_counts'], bins, alpha = 0.5, label = 'doublet')
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'hashtag'), 'n_counts'], bins, alpha = 0.5, label = 'singlet')
	ax.legend(loc='upper right')
	ax.set_xlabel("Number of UMIs")
	ax.set_ylabel("Number of nuclei")

	ax = axes[1, 1]
	ax.hist(adata.obs.loc[idx & np.isin(adata.obs['Cond2'], 'nohashtag'), 'n_counts'], bins, label = 'nohash')
	ax.legend(loc='upper right')
	ax.set_xlabel("Number of UMIs")
	ax.set_ylabel("Number of nuclei")

	fig.suptitle(nuclei_type)

	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()

def plot_doublet_hists(adata, out_file):
	with PdfPages(out_file) as pdf:
		for nuclei_type in adata.obs['coarse_annotation'].cat.categories:
			plot_doublet_hist(adata, nuclei_type, pdf, format = "pdf")

