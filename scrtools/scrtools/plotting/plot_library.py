import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted



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

