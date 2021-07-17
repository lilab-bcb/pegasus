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
  pegasus plot [options] [--restriction <restriction>...] [--palette <palette>...] <plot_type> <input_file> <output_file>
  pegasus plot -h

Arguments:
  plot_type              Plot type, either 'scatter' for scatter plots, 'compo' for composition plots, or 'wordcloud' for word cloud plots.
  input_file             Single cell data in Zarr or H5ad format.
  output_file            Output image file.

Options:
  --dpi <dpi>                        DPI value for the figure. [default: 500]
  --restriction <restriction>...     Set restriction if you only want to plot a subset of data. Multiple <restriction> strings are allowed. There are two types of <restriction>s: global restriction and attribute-specific restriction. Global restriction appiles to all attributes in ``attrs`` and takes the format of 'key:value,value...', or 'key:~value,value...'. This restriction selects cells with the ``data.obs[key]`` values belong to 'value,value...' (or not belong to if '~' shows). Attribute-specific restriction takes the format of 'attr:key:value,value...', or 'attr:key:~value,value...'. It only applies to one attribute 'attr'. If 'attr' and 'key' are the same, one can use '.' to replace 'key' (e.g. ``cluster_labels:.:value1,value2``).Each <restriction> takes the format of 'attr:value,value', or 'attr:~value,value..." which means excluding values. This option is used in both 'composition' and 'scatter'.
  --attributes <attrs>               <attrs> is a comma-separated list of attributes to color the basis. This option is only used in 'scatter'.
  --basis <basis>                    Basis for 2D plotting, chosen from 'umap', 'tsne', 'fitsne', 'fle', 'net_umap', 'net_tsne', 'net_fitsne' or 'net_fle'. [default: umap]
  --alpha <alpha>                    Point transparent parameter. Can be a single value or a list of values separated by comma used for each attribute in <attrs>.
  --legend-loc <str>                 Legend location, can be either "right margin" or "on data". If a list is provided, set 'legend_loc' for each attribute in 'attrs' separately. [default: right margin]
  --palette <str>                    Used for setting colors for every categories in categorical attributes. Multiple <palette> strings are allowed. Each string takes the format of 'attr:color1,color2,...,colorn'. 'attr' is the categorical attribute and 'color1' - 'colorn' are the colors for each category in 'attr' (e.g. 'cluster_labels:black,blue,red,...,yellow'). If there is only one categorical attribute in 'attrs', ``palletes`` can be set as a single string and the 'attr' keyword can be omitted (e.g. "blue,yellow,red").
  --show-background                  Show points that are not selected as gray.
  --nrows <nrows>                    Number of rows in the figure. If not set, pegasus will figure it out automatically.
  --ncols <ncols>                    Number of columns in the figure. If not set, pegasus will figure it out automatically.
  --panel-size <sizes>               Panel size in inches, w x h, separated by comma. Note that margins are not counted in the sizes. For composition, default is (6, 4). For scatter plots, default is (4, 4).
  --left <left>                      Figure's left margin in fraction with respect to panel width.
  --bottom <bottom>                  Figure's bottom margin in fraction with respect to panel height.
  --wspace <wspace>                  Horizontal space between panels in fraction with respect to panel width.
  --hspace <hspace>                  Vertical space between panels in fraction with respect to panel height.
  --groupby <attr>                   Use <attr> to categorize the cells for the composition plot, e.g. cell type.
  --condition <attr>                 Use <attr> to calculate frequency within each category defined by '--groupby' for the composition plot, e.g. donor.
  --style <style>                    Composition plot styles. Can be either 'frequency' or 'normalized'. [default: normalized]
  --factor <factor>                  Factor index (column index in data.uns['W']) to be used to generate word cloud plot.
  --max-words <max_words>            Maximum number of genes to show in the image. [default: 20]

  -h, --help                         Print out help information.

Examples:
  pegasus plot scatter --basis tsne --attributes louvain_labels,Donor example.h5ad scatter.pdf
  pegasus plot compo --groupby louvain_labels --condition Donor example.zarr.zip compo.pdf
  pegasus plot wordcloud --factor 0 example.zarr.zip word_cloud_0.pdf
    """
    def execute(self):
        kwargs = {
            "restrictions": self.args["--restriction"],
            "attrs": self.convert_to_list(self.args["--attributes"]),
            "basis": self.args["--basis"],
            "alpha": self.convert_to_list(self.args["--alpha"], converter=float),
            "legend_loc": self.convert_to_list(self.args["--legend-loc"]),
            "palettes" : self.args["--palette"],
            "show_background": self.args["--show-background"],
            "nrows": self.convert_to_int(self.args["--nrows"]),
            "ncols": self.convert_to_int(self.args["--ncols"]),
            "panel_size": self.convert_to_list(self.args["--panel-size"], converter=float),
            "left": self.convert_to_float(self.args["--left"]),
            "bottom": self.convert_to_float(self.args["--bottom"]),
            "wspace": self.convert_to_float(self.args["--wspace"]),
            "hspace": self.convert_to_float(self.args["--hspace"]),
            "groupby": self.args["--groupby"],
            "condition": self.args["--condition"],
            "style": self.args["--style"],
            "factor": int(self.args["--factor"]) if self.args["--factor"] is not None else self.args["--factor"],
            "max_words": int(self.args["--max-words"]),
            "return_fig": True,
            "dpi": int(self.args["--dpi"]),
        }

        for key in ["nrows", "ncols", "panel_size", "left", "bottom", "wspace", "hspace"]:
            if kwargs[key] is None:
                del kwargs[key]

        if self.args["<plot_type>"] == "scatter" and kwargs["attrs"] is None:
            raise KeyError("--attributes must be provided for scatter plots!")
        if self.args["<plot_type>"] == "compo" and (kwargs["groupby"] is None or kwargs["condition"] is None):
            raise KeyError("--groupby and --condition must be provided for composition plots!")
        if self.args["<plot_type>"] == "wordcloud" and kwargs["factor"] is None:
            raise KeyError("--factor must be provided for word cloud plots!")

        plot_type2keyword = {"scatter": "scatter", "compo" : "compo_plot", "wordcloud": "wordcloud"}

        data = read_input(self.args["<input_file>"])
        fig = getattr(pegasus.plotting, plot_type2keyword[self.args["<plot_type>"]])(data, **kwargs)

        output_file = self.args["<output_file>"]
        fig.savefig(output_file)
        logger.info(f"{output_file} is generated.")
