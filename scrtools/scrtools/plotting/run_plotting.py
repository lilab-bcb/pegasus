import scanpy.api as sc

from . import plot_library

def make_plots(input_file, plot_type, output_file, **kwargs):
	adata = sc.read(input_file)
	if plot_type == 'composition':
		plot_library.plot_composition(adata.obs['louvain_labels'], adata.obs[kwargs['attr']], output_file, style = kwargs['style'], stacked = kwargs['stacked'], logy = kwargs['logy'], sizes = kwargs['sizes'], rmove = kwargs['rmove'], wshrink = kwargs['wshrink'])
	elif plot_type == 'diffmap':
		assert len(kwargs['attrs']) > 0
		sc.pl.diffmap(adata, color = kwargs['attrs'], projection = kwargs['dim'], save = output_file)
