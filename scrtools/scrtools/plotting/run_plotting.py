import pandas as pd
from matplotlib import pyplot as pl

from ..tools import read_input
from .plot_utils import transform_basis
from . import plot_library, iplot_library


pop_list = {
	'composition' : {'basis', 'attrs', 'apply_to_all', 'group', 'genes', 'gene', 'nrows', 'ncols', 'alpha', 'legend_fontsize', 'use_raw', 'showzscore', 'title', 'showall'},
	'scatter' : {'cluster', 'attr', 'group', 'genes', 'gene', 'style', 'stacked', 'logy', 'use_raw', 'showzscore', 'title', 'showall'},
	'scatter_groups' : {'attr', 'attrs', 'apply_to_all', 'genes', 'gene', 'style', 'stacked', 'logy', 'use_raw', 'showzscore', 'title'},
	'scatter_genes' : {'cluster', 'attr', 'attrs', 'restrictions', 'apply_to_all', 'group', 'gene', 'style', 'stacked', 'logy', 'legend_fontsize', 'showzscore', 'title', 'showall'},
	'scatter_gene_groups' : {'cluster', 'attr', 'attrs', 'restrictions', 'apply_to_all', 'genes', 'style', 'stacked', 'logy', 'legend_fontsize', 'showzscore', 'title', 'showall'},
	'heatmap' : {'attr', 'basis', 'attrs', 'restrictions', 'apply_to_all', 'group', 'gene', 'style', 'stacked', 'logy', 'nrows', 'ncols', 'subplot_size', 'left', 'bottom', 'wspace', 'hspace', 'alpha', 'legend_fontsize', 'showall'}
}

def make_static_plots(input_file, plot_type, output_file, dpi = 500, **kwargs):
	adata = read_input(input_file, mode = 'r')
	assert plot_type in pop_list
	pop_set = pop_list[plot_type].copy()
	for key, value in kwargs.items():
		if value is None:
			pop_set.add(key)
	for key in pop_set:
		kwargs.pop(key)
	fig = getattr(plot_library, 'plot_' + plot_type)(adata, **kwargs)
	fig.savefig(output_file, dpi = dpi)
	print(output_file + " is generated.")


def make_interactive_plots(input_file, plot_type, output_file, **kwargs):
	adata = read_input(input_file, mode = 'r')
	basis = transform_basis(plot_type)
	if plot_type == 'diffmap' or plot_type == 'diffmap_pca':
		df = pd.DataFrame(adata.obsm['X_{}'.format(plot_type)][:, 0:3], index = adata.obs.index, columns = [basis + i for i in ['1', '2', '3']])
		if kwargs['isgene']:
			coln = adata.var.index.get_loc(kwargs['attr'])
			df.insert(0, 'Annotation', adata.X[:, coln].toarray().ravel())
		else:
			df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['isreal']:
			iplot_library.scatter3d(df, output_file)
		else:
			iplot_library.scatter3d_real(df, output_file, kwargs['log10'])
	else:
		df = pd.DataFrame(adata.obsm['X_{}'.format(plot_type)], index = adata.obs.index, columns = [basis + i for i in ['1', '2']])
		if kwargs['isgene']:
			coln = adata.var.index.get_loc(kwargs['attr'])
			df.insert(0, 'Annotation', adata.X[:, coln].toarray().ravel())
		else:
			df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['isreal']:
			iplot_library.scatter(df, output_file)
		else:
			iplot_library.scatter_real(df, output_file, kwargs['log10'])
	print(output_file + " is generated.")

