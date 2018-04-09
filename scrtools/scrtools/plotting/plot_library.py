import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted

import cycler
import plotly.offline as py
import plotly.graph_objs as go


def plot_composition(clust_labels, variable, outfile, style = 'frequency', stacked = True, logy = False, sizes = (6, 4), rmove = 1.1, wshrink = 0.8):
	fig, ax = plt.subplots(nrows = 1, ncols = 1)
	fig.set_size_inches(sizes)
	df = pd.crosstab(clust_labels, variable)
	df = df.reindex(index = natsorted(df.index.values), columns = natsorted(df.columns.values))

	if style == 'frequency':
		df = df.div(df.sum(axis = 1), axis = 0) * 100.0
	elif style == 'normalized':
		df = df.div(df.sum(axis = 0), axis = 1) * 100.0

	if logy and not stacked:
		df_sum = df.sum(axis = 1)
		df_new = df.cumsum(axis = 1)
		df_new = 10 ** df_new.div(df_sum, axis = 0).mul(np.log10(df_sum), axis = 0)
		df = df_new.diff(axis = 1).fillna(value = df_new.iloc[:, 0:1], axis = 1)
		df.plot(kind = 'bar', stacked = False, legend = False, logy = True, ylim = (1.01, df_sum.max() * 1.7), ax = ax)
	else:
		df.plot(kind = 'bar', stacked = stacked, legend = False, logy = logy, ax = ax)

	ax.set_xlabel('Cluster ID')
	ax.set_ylabel('Percentage' if style != 'count' else 'Count')
	ax.set_title("AMI = {0:.4f}".format(adjusted_mutual_info_score(clust_labels, variable)))
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * wshrink, box.height])
	ax.legend(loc = 'center', bbox_to_anchor = (rmove, 0.5), ncol = 1)
	ax.grid(False)
	fig.savefig(outfile)
	plt.close()



vega_20_scanpy = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
    '#9467bd', '#8c564b', '#e377c2',  # '#7f7f7f' removed grey
    '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896',
    '#c5b0d5', '#c49c94', '#f7b6d2',  # '#c7c7c7' removed grey
    '#dbdb8d', '#9edae5',
    '#ad494a', '#8c6d31']  # manual additions

palettes = cycler.cycler(color = vega_20_scanpy)

def plot_diffmap(df, output_file):
	labels = natsorted(df['Annotation'].unique())
	data = []
	for i, (label, color) in enumerate(zip(labels, palettes)):
		tdf = df[df['Annotation'] == label]
		trace = dict(
			name = label,
			x = tdf['x'],
			y = tdf['y'],
			z = tdf['z'],
			type = 'scatter3d',
			mode = 'markers',
			marker = dict(size = 5, color = color['color'], line = dict(width=0))
			)
		data.append(trace)
	layout = dict(
		title = 'test',
		scene = dict(
			xaxis = dict(title = 'DC1'),
			yaxis = dict(title = 'DC2'),
			zaxis = dict(title = 'DC3')
			),
		margin = dict(l = 0, r = 0, b = 0, t = 0)
		)
	fig = dict(data=data, layout=layout)
	py.plot(fig, filename = output_file)



def plot_kde(clust_labels, count_vector, clust_id, outfile, doublet_rate = 0.031, count_type = 'gene', format = None):
	fig, ax = plt.subplots(nrows = 1, ncols = 1)
	sr = pd.Series(count_vector) if clust_id == "all" else pd.Series(count_vector[clust_labels == clust_id])
	sr.plot(kind = 'kde', ax = ax)
	ax.axvline(x = sr.quantile(1.0 - doublet_rate), ls = '--', c = 'red')
	ax.grid(False)
	ax.tick_params(axis = 'x', labelsize = 8)
	ax.tick_params(axis = 'y', labelsize = 5)
	ax.set_xlabel("Number of {}s".format(count_type), fontsize = 10)
	ax.set_ylabel("Density", fontsize = 10)
	ax.set_title("Kernel density estimation of cluster {}".format(clust_id), fontsize = 12)
	plt.tight_layout(rect = (-0.05, 0, 1, 1))
	# outfile.savefig()
	plt.savefig(outfile, format = format)
	plt.close()

def plot_kde_for_all(clust_labels, count_vector, outfile, doublet_rate = 0.031, count_type = 'gene'):
	with PdfPages(outfile) as pdf:
		plot_kde(clust_labels, count_vector, "all", pdf, doublet_rate = doublet_rate, count_type = count_type, format = "pdf")
		print("done all.")
		for clust_id in natsorted(np.unique(clust_labels)):
			plot_kde(clust_labels, count_vector, clust_id, pdf, doublet_rate = doublet_rate, count_type = count_type, format = "pdf")
			print("done cluster {}".format(clust_id))


def load_table(inp_file):
	df = pd.read_table(inp_file, sep = '\t', header = 0, index_col = 0)
	for i, cn in enumerate(df):
		if i >= 2:
			if str(df[cn][0])[-1] == '%':
				df[cn] = df[cn].map(lambda x: float(x[:-1]))
			else:
				df[cn] = df[cn].map(lambda x: int(str(x).replace(',', '')))
	df.insert(0, 'Channel', df.index.values)

	return df


def make_2x2_plot(df, out_pdf, x, ys = ['Number of Reads', 'Sequencing Saturation', 'GRCh38 Median Percent Mito', 'GRCh38 Number Cells with High Quality'], hue = None, swarmplot = False, dot_size = 3):
	mpl.rcParams.update({'xtick.labelsize' : 5,
						 'ytick.labelsize' : 5,
						 'axes.labelsize' : 5})

	fig, axes = plt.subplots(nrows = 2, ncols = 2)
	plt.tight_layout(pad = 2)

	for i, ax in enumerate(axes.flatten()):		
		sns.boxplot(x = x, y = ys[i], hue = hue, data = df, orient = "v", fliersize = 0, dodge = True, linewidth = 1, ax = ax)
		if not swarmplot:
			sns.stripplot(x = x, y = ys[i], hue = hue, data = df, orient = "v", jitter = True, dodge = True, size = dot_size, linewidth = 0.5, edgecolor = 'gray', ax = ax)
		else:
			sns.stripplot(x = x, y = ys[i], hue = hue, data = df, orient = "v", dodge = True, size = dot_size, linewidth = 0.5, edgecolor = 'gray', ax = ax)

		if hue != None:	
			handles, labels = ax.get_legend_handles_labels()
			ax.legend(handles[:2], labels[:2])

		ax.set_xlabel("")
		a_list = []
		for tick in ax.get_xmajorticklabels():
			# a_list.append('_'.join(tick.get_text().split('_')[1:]))
			a_list.append('\n'.join(tick.get_text().split('_')))
		ax.set_xticklabels(a_list, rotation = 'vertical')
		# ax.tick_params(labelsize = 5)
		# ax.set_ylabel(ys[i], fontsize = 3)
		
	fig.savefig(out_pdf)
	plt.close()



def make_2x2_scatter_plot(dfHiSeq, dfNovaSeq, out_pdf, fields = ['Number of Reads', 'Sequencing Saturation', 'Reads Mapped Confidently to Transcriptome', 'Number Cells with High Quality']):
	mpl.rcParams.update({'font.size': 5})	

	fig, axes = plt.subplots(nrows = 2, ncols = 2)
	plt.tight_layout(pad = 5)

	for i, ax in enumerate(axes.flatten()):
		df = pd.DataFrame(data = {'HiSeq' : dfHiSeq[fields[i]].values, 'NovaSeq' : dfNovaSeq[fields[i]].values})
		df.plot.scatter(x = 'HiSeq', y = 'NovaSeq', ax = ax)
		lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
		ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
		ax.set_title(fields[i])
		
	fig.savefig(out_pdf)
	plt.close()


def load_qscore_table(inp_file, platform, ylab):
	df = pd.read_table(inp_file, sep = '\t', header = 0)
	df.insert(0, 'Platform', platform)
	df = df.melt(id_vars = ['Platform', 'Lane'], var_name = 'Position', value_name = ylab)
	if ylab != 'Qmean':
		df[ylab] = df[ylab].map(lambda x: x * 100.0)
	return df

def make_qscore_plot(path, rtype, out_pdf, sizes = (6, 4), platforms = ['NextSeq', 'HiSeq', 'NovaSeq'], plots = ['Qmean', 'Q30', 'Q20'], fontsize = 'small'):
	nr = len(plots)
	fig, axes = plt.subplots(nrows = nr, ncols = 1)
	fig.set_size_inches(sizes)
	
	if not isinstance(axes, np.ndarray):
		axes = np.array([axes])

	for i, (ax, ptype) in enumerate(zip(axes.flatten(), plots)):
		mylist = []
		for platform in platforms:
			mylist.append(load_qscore_table("{path}/{platform}_QC.{rtype}.{ptype}.txt".format(path = path, platform = platform, rtype = rtype, ptype = ptype), platform = platform, ylab = ptype))
		df = pd.concat(mylist)
		df['Position'] = df['Position'].map(int)
		sns.pointplot(x = 'Position', y = ptype, hue = 'Platform', data = df, capsize = 0.2, ax = ax, legend = False)
		ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.1), fancybox = True, ncol = len(platforms), fontsize = fontsize)

	fig.savefig(out_pdf)
	plt.close()
