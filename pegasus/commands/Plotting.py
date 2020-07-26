import os
from .Base import Base

import logging
logger = logging.getLogger(__name__)

from pegasusio import read_input
import pegasus.plotting



class Plotting(Base):
    """
Generate cluster composition plots.

Usage:
  pegasus plot [options] [--restriction <restriction>...] <plot_type> <input_file> <output_file>
  pegasus plot -h

Arguments:
  plot_type              Only 2D plots, chosen from 'composition', 'scatter'.
  input_file             Single cell data in Zarr or H5ad format.
  output_file            Output image file.

Options:
  --dpi <dpi>                        DPI value for the figure. [default: 500]
  --restriction <restriction>...     Set restriction if you only want to plot a subset of data. Multiple <restriction> strings are allowed. Each <restriction> takes the format of 'attr:value,value', or 'attr:~value,value..." which means excluding values. This option is used in 'composition' and 'scatter'.

  --attributes <attrs>               <attrs> is a comma-separated list of attributes to color the basis. This option is only used in 'scatter'.
  --basis <basis>                    Basis for 2D plotting, chosen from 'umap', 'tsne', 'fitsne', 'fle', 'net_umap', 'net_tsne', 'net_fitsne' or 'net_fle'. [default: umap]
  --alpha <alpha>                    Point transparent parameter. Can be a list of parameters separated by comma.
  --legend-loc <str>                 Legend location, can be either "right margin" or "on data" [default: right margin]
  --legend-fontsize <fontsize>       Legend font size.
  --apply-to-each-figure             Indicate that the <restriction> strings are not applied to all attributes but for specific attributes. The string's 'attr' value should math the attribute you want to restrict.
  --set-palettes <str>               A comma-separated list of colors for visualization.
  --show-background                  Show points that are not selected as gray.
  --nrows <nrows>                    Number of rows in the figure. If not set, pegasus will figure it out automatically.
  --ncols <ncols>                    Number of columns in the figure. If not set, pegasus will figure it out automatically.
  --subplot-size <sizes>             Sub-plot size in inches, w x h, separated by comma. Note that margins are not counted in the sizes. For composition, default is (6, 4). For scatter plots, default is (4, 4).
  --left <left>                      Figure's left margin in fraction with respect to subplot width.
  --bottom <bottom>                  Figure's bottom margin in fraction with respect to subplot height.
  --wspace <wspace>                  Horizontal space between subplots in fraction with respect to subplot width.
  --hspace <hspace>                  Vertical space between subplots in fraction with respect to subplot height.
  --do-not-show-all                  Do not show all components in group for scatter_groups.

  --xattr <attr>                     Use <attr> in x-axis for the composition plot, e.g. Donor.
  --yattr <attr>                     Use <attr> in y-axis for the composition plot, e.g. Cell type.
  --style <style>                    Composition plot styles. Can be either 'frequency', 'normalized'. [default: frequency]

  -h, --help                         Print out help information.

Examples:
  pegasus plot scatter --basis tsne --attributes louvain_labels,Individual Manton_BM.h5ad test.pdf
  pegasus plot composition --xattr Individual --yattr louvain_labels --style normalized Manton_BM.h5ad test.pdf
    """
    def execute(self):
        kwargs = {
            "restrictions": self.args["--restriction"],
            "attrs": self.convert_to_list(self.args["--attributes"]),
            "basis": self.args["--basis"],
            "alpha": self.convert_to_list(self.args["--alpha"], converter=float),
            "legend_loc": self.args["--legend-loc"],
            "legend_fontsize": self.convert_to_float(self.args["--legend-fontsize"]),
            "showall": not self.args["--do-not-show-all"],
            "palettes" : self.args["--set-palettes"],
            "show_background": self.args["--show-background"],
            "nrows": self.convert_to_int(self.args["--nrows"]),
            "ncols": self.convert_to_int(self.args["--ncols"]),
            "subplot_size": self.convert_to_list(self.args["--subplot-size"], converter=float),
            "left": self.convert_to_float(self.args["--left"]),
            "bottom": self.convert_to_float(self.args["--bottom"]),
            "wspace": self.convert_to_float(self.args["--wspace"]),
            "hspace": self.convert_to_float(self.args["--hspace"]),
            "xattr": self.args["--xattr"],
            "yattr": self.args["--yattr"],
            "style": self.args["--style"],
            "show": False,
            "dpi": int(self.args["--dpi"]),
        }

        for key in ["nrows", "ncols", "subplot_size", "left", "bottom", "wspace", "hspace"]:
            if kwargs[key] is None:
                del kwargs[key]

        plot_type2keyword = {"scatter": "scatter", "composition" : "compo_plot"}

        data = read_input(self.args["<input_file>"])
        fig = getattr(pegasus.plotting, plot_type2keyword[self.args["<plot_type>"]])(data, **kwargs)

        output_file = self.args["<output_file>"]
        fig.savefig(output_file)
        logger.info(f"{output_file} is generated.")
