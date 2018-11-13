import matplotlib as mpl
mpl.use("Agg")

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# plot_type: gene, count, mito
def plot_qc_violin(data, plot_type, out_file, xattr = 'Channel', hue = None, inner = None, dpi = 500, figsize = None, xlabel = None, xtick_font = None, xtick_rotation = False, no_legend = False):
	pt2attr = {'gene' : 'n_genes', 'count' : 'n_counts', 'mito' : 'percent_mito'}
	pt2ylab = {'gene' : 'Number of expressed genes', 'count' : 'Number of UMIs', 'mito' : 'Percentage of mitochondrial genes'}

	yattr = pt2attr[plot_type]

	split = True if hue is not None else False
	df = data.obs[[xattr, yattr]] if hue is None else data.obs[[xattr, yattr, hue]]

	sns.violinplot(x = xattr, y = yattr, hue = hue, data = df, inner = inner, split = split)
	ax = plt.gca()
	ax.grid(False)
	if xlabel is not None:
		ax.set_xlabel(xlabel)	
	ax.set_ylabel(pt2ylab[plot_type])
	for tick in ax.xaxis.get_major_ticks():
		if xtick_font is not None:
			tick.label.set_fontsize(xtick_font)
		if xtick_rotation: 
			tick.label.set_rotation('vertical')
	if no_legend:
		ax.legend_.remove()
	if figsize is not None:
		plt.gcf().set_size_inches(*figsize)
	plt.tight_layout()
	plt.savefig(out_file, dpi = dpi)
	plt.close()
