import numpy as np
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from natsort import natsorted
from scipy.stats import pearsonr

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



def plot_gene_violin(data, gene_list, out_file, dpi = 500, figsize = (6, 4), linewidth = None, quality = None):
	if quality is not None:
		data = data[data.obs['quality'] == quality]

	df = pd.DataFrame(data.X.toarray(), index = data.obs_names, columns = data.var_names)
	df['assignment'] = data.obs['demux_type'].astype(str)
	idx_singlet = np.isin(data.obs['demux_type'], 'singlet')
	df.loc[idx_singlet, 'assignment'] = data.obs.loc[idx_singlet, 'assignment'].astype(str)
	df['assignment'] = df['assignment'].astype('category')

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
			sns.violinplot(x = 'assignment', y = gene_list[gid], data = df, ax = ax, linewidth = linewidth)
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




def plot_human_vs_mouse(data, out_file, tp = 'counts', alpha = 1.0, dpi = 500, format = None, log = False):
	ax = plt.gca()

	if 'assignment' in data.obs:
		labels = data.obs['demux_type'].astype(str)
		labels[np.isin(data.obs['assignment'], ['1','2','3','4'])] = 'singlet_mouse'
		labels[np.isin(data.obs['assignment'], ['5','6','7','8'])] = 'singlet_human'
		labels = pd.Categorical(labels, categories = ['singlet_human', 'singlet_mouse', 'doublet', 'unknown'])
		colors = ['red', 'blue', 'green', 'orange']
	else:
		labels = ['nuclei'] * data.shape[0]
		labels = pd.Categorical(labels, categories = ['nuclei'])
		colors = ['blue']

	for k, cat in enumerate(labels.categories):
		idx = np.isin(labels, cat)
		ax.scatter(data.obs.loc[idx, 'mm10_n_' + tp],
				   data.obs.loc[idx, 'grch38_n_' + tp],
				   c = colors[k],
				   marker = '.',
				   alpha = alpha,
				   edgecolors = 'none',
				   label = cat,
				   rasterized = True)
	ax.grid(False)
	ax.legend()
	ax.set_xlabel("Mouse UMI")
	ax.set_ylabel("Human UMI")
	if log:
		ax.set_xscale('log')
		ax.set_yscale('log')
	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()


def plot_2attr_hvm(data, attrs, names, out_file, dpi = 500, format = None, log = True, separate = False, alpha = 0.5):	
	labels = data.obs['demux_type'].astype(str)
	labels[np.isin(data.obs['assignment'], ['1','2','3','4'])] = 'singlet_mouse'
	labels[np.isin(data.obs['assignment'], ['5','6','7','8'])] = 'singlet_human'
	labels = pd.Categorical(labels, categories = ['singlet_human', 'singlet_mouse', 'doublet', 'unknown'])
	colors = ['red', 'blue', 'green', 'orange']

	if separate:
		fig, axes = plt.subplots(nrows = 2, ncols = 2)
		axes = axes.flatten()
	else:
		ax = plt.gca()

	for k, cat in enumerate(labels.categories):
		idx = np.isin(labels, cat)

		if separate:
			ax = axes[k]

		x = data.obs.loc[idx, attrs[0]].values + 1.0
		y = data.obs.loc[idx, attrs[1]].values + 1.0

		ax.scatter(x,
				   y,
				   c = colors[k],
				   marker = '.',
				   alpha = alpha,
				   edgecolors = 'none',
				   label = cat,
				   rasterized = True)

		if separate:
			ax.grid(False)
			ax.set_title("{}, $\\log_{{10}}\\rho$ = {:.2f}".format(cat, pearsonr(np.log10(x), np.log10(y))[0]))
			if log:
				ax.set_xscale('log')
				ax.set_yscale('log')

	if separate:
		plt.tight_layout()
		plt.subplots_adjust(left = 0.1, bottom = 0.1)
		fig.text(0.5, 0.01, names[0], ha='center')
		fig.text(0.01, 0.5, names[1], va='center', rotation='vertical')
	else:
		ax.grid(False)
		ax.legend()
		x = data.obs[attrs[0]].values + 1.0
		y = data.obs[attrs[1]].values + 1.0
		ax.set_title("$\\log_{{10}}\\rho$ = {:.2f}".format(pearsonr(np.log10(x), np.log10(y))[0]))
		ax.set_xlabel(names[0])
		ax.set_ylabel(names[1])
		if log:
			ax.set_xscale('log')
			ax.set_yscale('log')


	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()


def plot_bar(heights, tick_labels, xlabel, ylabel, out_file, dpi = 500, figsize = None, format = None):
	plt.bar(x = np.linspace(0.5, heights.size - 0.5, heights.size), 
			height = heights, 
			tick_label = tick_labels)
	ax = plt.gca()
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if figsize is not None:
		fig = plt.gcf()
		fig.set_size_inches(*figsize)
	plt.tick_params(axis = 'x', labelsize = 7)
	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()

# attrs is a dict with name: attr format
def plot_violins(data, attrs, xlabel, ylabel, out_file, dpi = 500, format = None, log = True):
	dfs = []
	if isinstance(data, anndata.base.AnnData):
		for name, attr in attrs.items():
			dfs.append(pd.DataFrame({xlabel : name, ylabel : data.obs[attr].values}))
	else:
		for arr, name in zip(data, attrs):
			dfs.append(pd.DataFrame({xlabel : name, ylabel : arr}))
	df = pd.concat(dfs)

	if log:
		df[ylabel] = np.log10(df[ylabel])
		sns.violinplot(x = xlabel, y = ylabel, data = df)
		y_max = int(np.ceil(df[ylabel].max()))
		loc = list(range(y_max + 1))
		labels = [r'$10^{}$'.format(x) for x in loc]
		plt.yticks(loc, labels)
	else:
		sns.violinplot(x = xlabel, y = ylabel, data = df)

	plt.gca().grid(False)
	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()

def plot_demux_hist(data, out_file, plot_attr = 'n_counts', cat_attr = 'demux_type', dpi = 500, format = None):
	bins = np.logspace(np.log10(min(data.obs[plot_attr])), np.log10(max(data.obs[plot_attr])), 101)
	cat_vec = data.obs[cat_attr]
	ax = plt.gca()
	if cat_attr == 'demux_type':
		ax.hist(data.obs.loc[np.isin(cat_vec, 'singlet'), plot_attr], bins, alpha = 0.5, label = 'singlet')
		ax.hist(data.obs.loc[np.isin(cat_vec, 'doublet'), plot_attr], bins, alpha = 0.5, label = 'doublet')
		ax.hist(data.obs.loc[np.isin(cat_vec, 'unknown'), plot_attr], bins, alpha = 0.5, label = 'unknown')
	else:
		assert cat_attr == 'quality'
		ax.hist(data.obs.loc[np.isin(cat_vec, 'high'), plot_attr], bins, alpha = 0.5, label = 'high')
		ax.hist(data.obs.loc[np.isin(cat_vec, 'low'), plot_attr], bins, alpha = 0.5, label = 'low')
	ax.legend(loc='upper right')
	ax.set_xscale("log")
	ax.set_xlabel("Number of RNA UMIs (log10 scale)")
	ax.set_ylabel("Number of cellular barcodes")
	plt.savefig(out_file, dpi = dpi, format = format)
	plt.close()

def plot_heatmap(vec1, vec2, out_file, dpi = 500, format = None):
	df = pd.crosstab(vec1, vec2)
	df.columns.name = df.index.name = ''

	ax = plt.gca()
	ax.xaxis.tick_top()
	ax.set_xlabel('')
	ax.set_ylabel('')
	ax = sns.heatmap(df, annot = True, fmt = 'd', cmap = 'inferno', ax = ax)
	
	plt.tight_layout()
	plt.savefig(out_file, dpi = 500, format = format)
	plt.close()
