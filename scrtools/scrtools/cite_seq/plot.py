import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_barcode_hist(data, adt, out_file, dpi = 500):
	background = adt.X.sum(axis = 1).A1
	signal = adt[data.obs_names,].X.sum(axis = 1).A1
	bins = np.logspace(0, np.log10(max(background)), 100)
	plt.hist(background, bins, alpha = 0.5, label = 'background', log = True)
	plt.hist(signal, bins, alpha = 0.5, label = 'signal', log = True)
	plt.legend(loc='upper right')
	ax = plt.gca()
	ax.set_xscale("log")
	ax.set_xlabel("Number of UMIs (log10 scale)")
	ax.set_ylabel("Number of barcodes (log10 scale)")
	plt.savefig(out_file, dpi = dpi)
	plt.close()



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
	fig.savefig(out_file, dpi = dpi)

