import pandas as pd
import scanpy.api as sc
from scanpy.plotting.anndata import scatter
from matplotlib import pyplot as pl

from . import plot_library

def make_static_plots(input_file, plot_type, output_name, **kwargs):
	adata = sc.read(input_file)
	if plot_type == 'composition':
		output_file = output_name + '.' + kwargs['ext']
		plot_library.plot_composition(adata.obs['louvain_labels'], adata.obs[kwargs['attr']], output_file, style = kwargs['style'], stacked = kwargs['stacked'], logy = kwargs['logy'], sizes = kwargs['sizes'], rmove = kwargs['rmove'], wshrink = kwargs['wshrink'])
	elif plot_type == 'diffmap':
		assert len(kwargs['attrs']) > 0
		plot_diffusion_map(adata, kwargs['attrs'], output_name, ext = kwargs['ext'], 
			which = kwargs['angles'], projection = kwargs['dim'], legend_fontsize = kwargs['lfsize'], size = kwargs['ptsize'])

def make_interactive_plots(input_file, plot_type, output_file, **kwargs):
	adata = sc.read(input_file)
	if plot_type == 'diffmap':
		df = pd.DataFrame(adata.obsm['X_diffmap'][:, 0:3], index = adata.obs.index, columns = ['x', 'y', 'z'])
		df.insert(0, 'Annotation', adata.obs[kwargs['attr']])
		if not kwargs['real']:
			plot_library.plot_diffmap(df, output_file)
		else:
			plot_library.plot_diffmap_real(df, output_file, kwargs['log10'])

def plot_diffusion_map(adata, color, output_name, ext = 'png', which = 'all', projection = '3d', legend_loc = 'right margin', legend_fontsize = None, legend_fontweight = None, size = None, title = None):
	proj2d = ['1,2', '2,3', '1,3']
	proj3d = ['1,2,3', '2,3,1', '1,3,2']

	if which == 'all':
		which = '0,1,2'
	proj = proj2d if projection == '2d' else proj3d
	which = [int(x) for x in which.split(',')]
	for cid in which:
		components = proj[cid]
		axs = scatter(
			adata,
			basis='diffmap',
			color=color,
			components=components,
			projection=projection,
			legend_loc=legend_loc,
			legend_fontsize=legend_fontsize,
			legend_fontweight=legend_fontweight,
			size=size,
			title=title,
			show=False,
			save=False)
		output_file = "{0}_{1}.{2}".format(output_name, components, ext)
		pl.savefig(output_file)
		print("{0} is generated.".format(output_file))
