import matplotlib as mpl
mpl.use("Agg")

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# plot_type: gene, count, mito
def plot_qc_violin(data, plot_type, out_file, xattr = 'Channel', hue = None, inner = None, dpi = 500, figsize = None, xlabel = None, xtick_font = None, xtick_rotation = False, split = False, linewidth = None):
	pt2attr = {'gene' : 'n_genes', 'count' : 'n_counts', 'mito' : 'percent_mito'}
	pt2ylab = {'gene' : 'Number of expressed genes', 'count' : 'Number of UMIs', 'mito' : 'Percentage of mitochondrial UMIs'}

	yattr = pt2attr[plot_type]

	tmp_df = data if isinstance(data, pd.core.frame.DataFrame) else data.obs
	df = tmp_df[[xattr, yattr]] if hue is None else tmp_df[[xattr, yattr, hue]]

	sns.violinplot(x = xattr, y = yattr, hue = hue, data = df, inner = inner, split = split, linewidth = linewidth, cut = 0)
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
	ax.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = xtick_font)
	if figsize is not None:
		plt.gcf().set_size_inches(*figsize)
	plt.tight_layout()
	plt.savefig(out_file, dpi = dpi)
	plt.close()
