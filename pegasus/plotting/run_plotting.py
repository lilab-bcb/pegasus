import pandas as pd
from matplotlib import pyplot as pl

from pegasusio import read_input
from .plot_utils import transform_basis
from .plot_qc import plot_qc_violin
from . import plot_library




def make_static_plots(input_file, plot_type, output_file, dpi=500, **kwargs):
    adata = read_input(input_file, mode="r")

    if plot_type == "qc_violin":
        if kwargs["attr"] is None:
            plot_qc_violin(
                adata,
                kwargs["qc_type"],
                output_file,
                xattr=kwargs["cluster"],
                xlabel=kwargs["cluster"],
                xtick_font=kwargs["qc_xtick_font"],
                xtick_rotation=kwargs["qc_xtick_rotation"],
                figsize=kwargs["subplot_size"],
                linewidth=kwargs["qc_line_width"],
            )
        else:
            plot_qc_violin(
                adata,
                kwargs["qc_type"],
                output_file,
                xattr=kwargs["cluster"],
                hue=kwargs["attr"],
                xlabel=kwargs["cluster"],
                xtick_font=kwargs["qc_xtick_font"],
                xtick_rotation=kwargs["qc_xtick_rotation"],
                split=True,
                figsize=kwargs["subplot_size"],
                linewidth=kwargs["qc_line_width"],
            )
    else:
        fig = getattr(plot_library, "plot_" + plot_type)(adata, **kwargs)
        fig.savefig(output_file, dpi=dpi)

    print(output_file + " is generated.")

