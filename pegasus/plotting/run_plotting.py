import pandas as pd
from matplotlib import pyplot as pl

from pegasus.io import read_input
from .plot_utils import transform_basis
from .plot_qc import plot_qc_violin
from . import plot_library, iplot_library




def make_static_plots(input_file, plot_type, output_file, dpi=500, **kwargs):
    adata = read_input(input_file, h5ad_mode="r")

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
    adata.file.close()


def make_interactive_plots(input_file, plot_type, output_file, **kwargs):
    adata = read_input(input_file, h5ad_mode="r")
    basis = transform_basis(plot_type)
    if plot_type == "diffmap" or plot_type == "diffmap_pca":
        df = pd.DataFrame(
            adata.obsm["X_{}".format(plot_type)][:, 0:3],
            index=adata.obs.index,
            columns=[basis + i for i in ["1", "2", "3"]],
        )
        if kwargs["isgene"]:
            coln = adata.var.index.get_loc(kwargs["attr"])
            df.insert(0, "Annotation", adata.X[:, coln].toarray().ravel())
        else:
            df.insert(0, "Annotation", adata.obs[kwargs["attr"]])
        if not kwargs["isreal"]:
            iplot_library.scatter3d(df, output_file)
        else:
            iplot_library.scatter3d_real(df, output_file, kwargs["log10"])
    else:
        df = pd.DataFrame(
            adata.obsm["X_{}".format(plot_type)],
            index=adata.obs.index,
            columns=[basis + i for i in ["1", "2"]],
        )
        if kwargs["isgene"]:
            coln = adata.var.index.get_loc(kwargs["attr"])
            df.insert(0, "Annotation", adata.X[:, coln].toarray().ravel())
        else:
            df.insert(0, "Annotation", adata.obs[kwargs["attr"]])
        if not kwargs["isreal"]:
            iplot_library.scatter(df, output_file)
        else:
            iplot_library.scatter_real(df, output_file, kwargs["log10"])
    print(output_file + " is generated.")
    adata.file.close()
