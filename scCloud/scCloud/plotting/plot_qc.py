import matplotlib as mpl
mpl.use("Agg")

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# plot_type: gene, count, mito
def plot_qc_violin(data, plot_type, out_file, dpi = 500, figsize = None):
	pt2attr = {'gene' : 'n_genes', 'count' : 'n_counts', 'mito' : 'percent_mito'}
	pt2ylab = {'gene' : 'Number of expressed genes', 'count' : 'Number of counts', 'mito' : 'Percentage of mitochondrial genes'}

	attr = pt2attr[plot_type]
	sns.violinplot(x = 'Channel', y = attr, data = data.obs[['Channel', attr]])
	sns.stripplot(x = 'Channel', y = attr, data = data.obs[['Channel', attr]], orient = 'v', color = 'black')
	ax = plt.gca()
	ax.grid(False)
	ax.set_ylabel(pt2ylab[plot_type])
	if figsize is not None:
		plt.gcf().set_size_inches(*figsize)
	plt.tight_layout()
	plt.savefig(out_file, dpi = dpi)
	plt.close()
