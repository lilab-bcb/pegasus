import pandas as pd
import anndata
from matplotlib import pyplot as pl

from .plot_utils import transform_basis
from . import plot_library, iplot_library



pop_list = {
	'composition' : {'basis', 'attrs', 'group', 'genes', 'gene', 'nrows', 'ncols', 'alpha', 'legend_fontsize', 'use_raw', 'showzscore', 'title'},
	'scatter' : {'cluster', 'attr', 'group', 'genes', 'gene', 'style', 'stacked', 'logy', 'use_raw', 'showzscore', 'title'},
	'scatter_groups' : {'attr', 'attrs', 'genes', 'gene', 'style', 'stacked', 'logy', 'use_raw', 'showzscore', 'title'},
	'scatter_genes' : {'cluster', 'attr', 'attrs', 'group', 'gene', 'style', 'stacked', 'logy', 'legend_fontsize', 'showzscore', 'title'},
	'scatter_gene_groups' : {'cluster', 'attr', 'attrs', 'genes', 'style', 'stacked', 'logy', 'legend_fontsize', 'showzscore', 'title'},
	'heatmap' : {'attr', 'basis', 'attrs', 'group', 'gene', 'style', 'stacked', 'logy', 'nrows', 'ncols', 'subplot_size', 'left', 'bottom', 'wspace', 'hspace', 'alpha', 'legend_fontsize'}
}

def make_static_plots(input_file, plot_type, output_file, **kwargs):
	adata = anndata.read_h5ad(input_file)
	print("Input file is loaded.")	
	assert plot_type in pop_list
	pop_set = pop_list[plot_type].copy()
	for key, value in kwargs.items():
		if value is None:
			pop_set.add(key)
	for key in pop_set:
		kwargs.pop(key)
	print(kwargs)
	fig = getattr(plot_library, 'plot_' + plot_type)(adata, **kwargs)
	fig.savefig(output_file)


def make_interactive_plots(input_file, plot_type, output_file, **kwargs):
	adata = anndata.read_h5ad(input_file)
	print("Input file is loaded.")	
	basis = transform_basis(plot_type)
	if plot_type == 'diffmap':
		df = pd.DataFrame(adata.obsm['X_diffmap'][:, 0:3], index = adata.obs.index, columns = [basis + i for i in ['1', '2', '3']])
		if kwargs['isgene']:
			df.insert(0, 'Annotation', adata[:, kwargs['attr']].X)
		else:
			df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['isreal']:
			iplot_library.scatter3d(df, output_file)
		else:
			iplot_library.scatter3d_real(df, output_file, kwargs['log10'])
	else:
		df = pd.DataFrame(adata.obsm['X_{}'.format(plot_type)], index = adata.obs.index, columns = [basis + i for i in ['1', '2']])
		if kwargs['isgene']:
			df.insert(0, 'Annotation', adata[:, kwargs['attr']].X)
		else:
			df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['isreal']:
			iplot_library.scatter(df, output_file)
		else:
			iplot_library.scatter_real(df, output_file, kwargs['log10'])
