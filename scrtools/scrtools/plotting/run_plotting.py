import pandas as pd
import anndata
from matplotlib import pyplot as pl

from . import plot_library

def make_static_plots(input_file, plot_type, output_name, **kwargs):
	adata = anndata.read_h5ad(input_file)
	if plot_type == 'composition':
		output_file = output_name + '.' + kwargs['ext']
		plot_library.plot_composition(adata.obs['louvain_labels'], adata.obs[kwargs['attr']], output_file, style = kwargs['style'], stacked = kwargs['stacked'], logy = kwargs['logy'], sizes = kwargs['sizes'], rmove = kwargs['rmove'], wshrink = kwargs['wshrink'])
	elif plot_type == 'diffmap':
		assert len(kwargs['attrs']) > 0
		plot_diffusion_map(adata, kwargs['attrs'], output_name, ext = kwargs['ext'], 
			which = kwargs['angles'], projection = kwargs['dim'], legend_fontsize = kwargs['lfsize'], size = kwargs['ptsize'])

def make_interactive_plots(input_file, plot_type, output_file, **kwargs):
	adata = anndata.read_h5ad(input_file)
	if plot_type == 'diffmap':
		df = pd.DataFrame(adata.obsm['X_diffmap'][:, 0:3], index = adata.obs.index, columns = ['DC1', 'DC2', 'DC3'])
		df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['real']:
			iplot_library.scatter3d(df, output_file)
		else:
			iplot_library.scatter3d_real(df, output_file, kwargs['log10'])
	else:
		df = pd.DataFrame(adata.obsm['X_{}'.format(plot_type)], index = adata.obs.index, columns = ['x', 'y'])
		df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['real']:
			iplot_library.scatter(df, output_file)
		else:
			iplot_library.scatter_real(df, output_file, kwargs['log10'])
