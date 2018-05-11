import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import zscore
from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted



def composition_plot(clust_labels, variable, style = 'frequency', stacked = True, logy = False, sizes = (6, 4), rmove = 1.1, wshrink = 0.8):
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

	return fig



### Extended color list for large datasets with many clusters
flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
colors_extra = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
mycols = sns.color_palette() + sns.color_palette("Set2", 10) + sns.color_palette(flatui) + sns.xkcd_palette(colors_extra) + sns.color_palette("hls", 8)
mypal = sorted(set(mycols), key=lambda c: mycols.index(c))
#sns.palplot(mypal) # disply colors


### Sample usage:
###    fig = scatter_plot(data, basis='tsne', components=[1,2], figsize=(15,15), title="Markers")
###    fig = scatter_plot(data_g, 'tsne', [1,2], genes=['PPBP','PF4'], nrows=1, ncols=2, figsize=(8,3))
def scatter_plot(data, basis, genes, components=[1,2], nrows=1, ncols=1, use_raw=False, title="", figsize=(15,15)):
	components = np.array(components).astype(int) - 1
	df = pd.DataFrame(data.obsm['X_' + basis][:, components], columns=[basis + str(c) for c in components])
	ngenes = len(genes)

	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=1, hspace=None)
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols, frameon=True, sharex=True, sharey=True, figsize=figsize)
	for i in range(nrows):
		for j in range(ncols):
			if ngenes <= (i*ncols + j):
				break

			if ncols == 1:
				ax = axes[i]
			elif nrows==1:
				ax = axes[j]
			else:
				ax = axes[i,j]

			marker = genes[i*ncols + j]
			ax.scatter(df[basis + str(components[0])].values,
					   df[basis + str(components[1])].values,
					   s=1,
					   c=data.raw[:,marker].X if use_raw else data[:,marker].X,
					   cmap='viridis',
					   edgecolor='none')
			ax.grid(False)
			### add the colorbar (this is what pandas does ref: https://github.com/pandas-dev/pandas/blob/3e160d21d2ac8dd56c008387a3d52d9f6a4d0f4d/pandas/plotting/_core.py#L873)
			img = ax.collections[0]
			kws = dict(ax=ax)
			ax.get_figure().colorbar(img, **kws)
			### set title and axis titles
			ax.set_title(marker)
			if i==(nrows-1): ax.set_xlabel(basis + str(components[0]))
			if j==0: ax.set_ylabel(basis + str(components[1]))

	return fig



def scatter_plot_cat(data, basis, attrs, components = [1, 2], nrows = 1, ncols = 1, figsize = (15, 15), rmove = 1.1, wshrink = 0.8):
	components = np.array(components).astype(int) - 1
	df_basis = pd.DataFrame(data.obsm['X_' + basis][:, components], columns = [basis + str(c) for c in components])
	nattrs = len(attrs)
	
	colorlist = np.array(mypal)

	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=1, hspace=None)
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols, frameon=True, sharex=True, sharey=True, figsize=figsize)
	for i in range(nrows):
		for j in range(ncols):
			if nattrs <= (i * ncols + j):
				break

			if ncols == 1:
				ax = axes[i]
			elif nrows == 1:
				ax = axes[j]
			else:
				ax = axes[i, j]

			attr = attrs[i * ncols + j]
			labels = data.obs[attr].astype('category')
			size = labels.cat.categories.size
			cat2color = dict(zip(labels.cat.categories.values, colorlist[range(size),:]))
			colors = labels.apply(lambda x: cat2color[x])
			colors.name = 'color'
			df = pd.concat([df_basis, labels, colors])
			df.plot(kind = 'scatter', x = 0, y = 1, c = 'color', label = attr, ax = ax)
			
			ax.grid(False)	
			ax.set_title(attr)
			box = ax.get_position()
			ax.set_position([box.x0, box.y0, box.width * wshrink, box.height])
			ax.legend(loc = 'center', bbox_to_anchor = (rmove, 0.5))
			if i == nrows - 1:
				ax.set_xlabel(basis + str(components[0]))
			if j == 0 :
				ax.set_ylabel(basis + str(components[1]))

	return fig

### Sample usage:
###     sns.set(font_scale=0.35)
###     cg = heatmap_plot(data_g, pd.DataFrame(adata_g.uns['rank_genes_groups_gene_names']), 20, use_raw=True, title="markers", figsize=(15, 7))
###     cg.savefig("heatmap.png", bbox_inches='tight', dpi=600)
def heatmap_plot(adata, ranked_genes, count, use_raw=False, showzscore=False, title="", **kwargs):
	if(isinstance(ranked_genes, np.recarray)): # e.g. adata.uns['rank_genes_groups_gene_names']
		ranked_genes = pd.DataFrame(ranked_genes)
	# e.g. pd.DataFrame(adata.uns['rank_genes_groups_gene_names'])
	genes = np.transpose(ranked_genes.iloc[0:count].values).flatten()
	genes = np.array(genes, dtype='str')
	genes = list(genes)
	genes = sorted(set(genes), key=lambda g: genes.index(g)) # dedup but keep order - does not actually matter since we cluster

	if use_raw:
		df = pd.DataFrame(adata.raw[:,genes].X.toarray(), index=adata.obs.index, columns=genes)
	else:
		df = pd.DataFrame(adata.X[:genes], index=adata.obs.index, columns=genes)

	if showzscore:
		df = df.apply(zscore, axis=0)

	df['Group'] = adata.smp['louvain_groups']
	ordered = df.sort_values('Group')

	unique = df['Group'].unique()
	palette = dict(zip(unique, mypal[0:len(unique)]))
	row_colors = df['Group'].map(palette)
	ordered = ordered.drop('Group', axis=1)

	cg = sns.clustermap(data=ordered,
				   cmap='RdBu_r',
				   row_colors=row_colors,
				   row_cluster=False, col_cluster=True,
				   linewidths=0,
				   yticklabels=[], xticklabels=genes, **kwargs)
	cg.ax_heatmap.set_ylabel('')
	# move the colorbar
	cg.ax_row_dendrogram.set_visible(False)
	dendro_box = cg.ax_row_dendrogram.get_position()
	dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
	dendro_box.x1 = dendro_box.x0 + 0.02
	cg.cax.set_position(dendro_box)
	cg.cax.yaxis.set_ticks_position("left")
	cg.cax.tick_params(labelsize=10)
	# draw a legend for the cluster groups
	cg.ax_col_dendrogram.clear()
	for label in unique:
		cg.ax_col_dendrogram.bar(0, 0, color=palette[label],label=label, linewidth=0)
	cg.ax_col_dendrogram.legend(loc="center", ncol=15, fontsize=10)
	cg.ax_col_dendrogram.grid(False)
	cg.ax_col_dendrogram.set_xticks([])
	cg.ax_col_dendrogram.set_yticks([])

	return cg
