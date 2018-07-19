import os
import time
from .Base import Base
from ..plotting import make_static_plots

class Plotting(Base):
    """
Generate cluster composition plots.

Usage:
  scrtools plot [options] [--restriction <restriction>...] <plot_type> <input_h5ad_file> <output_file>
  scrtools plot -h

Arguments:
  plot_type              Only 2D plots, chosen from 'composition', 'scatter', 'scatter_groups', 'scatter_genes', 'scatter_gene_groups', and 'heatmap'.
  input_h5ad_file        Single cell data with clustering done by Scanpy in h5ad file format.
  output_file            Output image file.

Options:
  --dpi <dpi>                        DPI value for the figure. [default: 500]

  --cluster-labels <attr>            Use <attr> as cluster labels. This option is used in 'composition', 'scatter_groups', and 'heatmap'.
  --attribute <attr>                 Plot <attr> against cluster labels. This option is only used in 'composition'.
  --basis <basis>                    Basis for 2D plotting, chosen from 'tsne', 'fitsne', 'umap', 'pca', 'rpca', 'fle', or 'diffmap_pca'. This option is used in 'scatter', 'scatter_groups', 'scatter_genes', and 'scatter_gene_groups'. [default: tsne]
  --attributes <attrs>               <attrs> is a comma-separated list of attributes to color the basis. This option is only used in 'scatter'.
  --restriction <restriction>...     Multiple <restriction> strings for different attributes. Each <restriction> takes the format of 'attr:value,value'. Only used for scatter. 
  --group <attr>                     <attr> is used to make group plots. In group plots, the first one contains all components in the group and the following plots show each component separately. This option is iused in 'scatter_groups' and 'scatter_gene_groups'. If <attr> is a semi-colon-separated string, parse the string as groups.
  --genes <genes>                    <genes> is a comma-separated list of gene names to visualize. This option is used in 'scatter_genes' and 'heatmap'.
  --gene <gene>                      Visualize <gene> in group plots. This option is only used in 'scatter_gene_groups'.

  --style <style>                    Composition plot styles. Can be either 'frequency', 'count', or 'normalized'. [default: frequency]
  --not-stacked                      Do not stack bars in composition plot.
  --log-y                            Plot y axis in log10 scale for composition plot.

  --nrows <nrows>                    Number of rows in the figure. If not set, scrtools will figure it out automatically.
  --ncols <ncols>                    Number of columns in the figure. If not set, scrtools will figure it out automatically.
  --subplot-size <sizes>             Sub-plot size in inches, w x h, separated by comma. Note that margins are not counted in the sizes. For composition, default is (6, 4). For scatter plots, default is (4, 4). 
  --left <left>                      Figure's left margin in fraction with respect to subplot width.
  --bottom <bottom>                  Figure's bottom margin in fraction with respect to subplot height.
  --wspace <wspace>                  Horizontal space between subplots in fraction with respect to subplot width.
  --hspace <hspace>                  Vertical space between subplots in fraction with respect to subplot height.
  --alpha <alpha>                    Point transparent parameter.
  --legend-fontsize <fontsize>       Legend font size.
  --use-raw                          Use anndata stored raw expression matrix. Only used by 'scatter_genes' and 'scatter_gene_groups'.

  --do-not-show-all                  Do not show all components in group for scatter_groups.

  --show-zscore                      If show zscore in heatmap.
  --heatmap-title <title>            Title for heatmap.

  -h, --help                         Print out help information.

Examples:
  scrtools plot composition --cluster-labels louvain_labels --attribute Individual --style normalized --not-stacked Manton_BM.h5ad test.png
  scrtools plot scatter --basis tsne --attributes louvain_labels,Individual Manton_BM.h5ad test.png
  scrtools plot scatter_groups --cluster-labels louvain_labels --group Individual Manton_BM.h5ad test.png
  scrtools plot scatter_genes --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP Manton_BM.h5ad test.png
  scrtools plot scatter_gene_groups --gene CD8A --group Individual Manton_BM.h5ad test.png
  scrtools plot heatmap --cluster-labels louvain_labels --genes CD8A,CD4,CD3G,MS4A1,NCAM1,CD14,ITGAX,IL3RA,CD38,CD34,PPBP --heatmap-title 'markers' Manton_BM.h5ad test.png
    """
    def execute(self):
        kwargs = {
            'cluster' : self.args['--cluster-labels'],
            'attr' : self.args['--attribute'],
            'restrictions' : self.args['--restriction'],
            'basis' : self.args['--basis'],
            'attrs' : self.split_string(self.args['--attributes']),
            'group' : self.args['--group'],
            'genes' : self.split_string(self.args['--genes']),
            'gene' : self.args['--gene'],
            'style' : self.args['--style'],
            'stacked' : not self.args['--not-stacked'],
            'logy' : self.args['--log-y'],
            'nrows' : int(self.args['--nrows']) if self.args['--nrows'] is not None else None,
            'ncols' : int(self.args['--ncols']) if self.args['--ncols'] is not None else None,
            'subplot_size' : [float(x) for x in self.args['--subplot-size'].split(',')] if self.args['--subplot-size'] is not None else None,
            'left' : float(self.args['--left']) if self.args['--left'] is not None else None,
            'bottom' : float(self.args['--bottom']) if self.args['--bottom'] is not None else None,
            'wspace' : float(self.args['--wspace']) if self.args['--wspace'] is not None else None,
            'hspace' : float(self.args['--hspace']) if self.args['--hspace'] is not None else None,
            'alpha' : float(self.args['--alpha']) if self.args['--alpha'] is not None else None,
            'legend_fontsize' : float(self.args['--legend-fontsize']) if self.args['--legend-fontsize'] is not None else None,
            'use_raw' : self.args['--use-raw'],
            'showzscore' : self.args['--show-zscore'],
            'title' : self.args['--heatmap-title'],
            'showall' : not self.args['--do-not-show-all']
        }
        
        make_static_plots(self.args['<input_h5ad_file>'], self.args['<plot_type>'], self.args['<output_file>'], dpi = int(self.args['--dpi']), **kwargs)
        time.sleep(1) # wait 1s to release the backed h5ad file
