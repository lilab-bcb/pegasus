import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.sparse import issparse
from pandas.api.types import is_numeric_dtype, is_list_like, is_categorical

from scipy.stats import zscore
from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted
from matplotlib import rcParams


import anndata
from pegasusio import UnimodalData, MultimodalData

from typing import List, Tuple, Union, Optional

import logging
logger = logging.getLogger(__name__)

from .plot_utils import _transform_basis, _get_nrows_and_ncols, _get_marker_size, _get_subplot_layouts, _get_legend_ncol, _get_palettes, RestrictionParser



def scatter(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attrs: Union[str, List[str]],
    basis: Optional[str] = "umap",
    matkey: Optional[str] = None,
    alpha: Optional[Union[float, List[float]]] = 1.0,
    legend_loc: Optional[str] = "right margin",
    legend_fontsize: Optional[int] = None,
    apply_to_all: Optional[bool] = True,
    restrictions: Optional[List[str]] = None,
    palettes: Optional[str] = None,
    show_background: Optional[bool] = False,
    legend_ncol: Optional[str] = None,
    cmap: Optional[str] = "YlOrRd",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    subplot_size: Optional[Tuple[float, float]] = (4, 4),
    left: Optional[float] = 0.2,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.4,
    hspace: Optional[float] = 0.15,
    show: Optional[bool] = True,
    **others,
) -> plt.Figure:
    """
    attrs: can be either in data.obs or a gene name
    legend_loc can be either "right margin" or  "on data"
    Example:
    pg.scatter(data, 'umap', 'leiden_labels')
    fig = pg.scatter(data, 'tsne', ['louvain_labels', 'hdbscan_labels_soft'], nrows = 1, ncols = 2, alpha = 0.5, show = False)
    """
    if not is_list_like(attrs):
        attrs = [attrs]
    nattrs = len(attrs)
    if not is_list_like(alpha):
        alpha = [alpha] * nattrs
    if matkey is not None:
        data.select_matrix(matkey)

    x = data.obsm[f"X_{basis}"][:, 0]
    y = data.obsm[f"X_{basis}"][:, 1]

    basis = _transform_basis(basis)
    marker_size = _get_marker_size(x.size)
    nrows, ncols = _get_nrows_and_ncols(nattrs, nrows, ncols)
    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=ncols, subplot_size=subplot_size, left=left, bottom=bottom, wspace=wspace, hspace=hspace, squeeze=False)

    if legend_loc is None:
        legend_loc = rcParams["legend.loc"]
    if legend_fontsize is None:
        legend_fontsize = rcParams["legend.fontsize"]

    restr_obj = RestrictionParser(restrictions)
    unsel = restr_obj.get_unsatisfied(data, apply_to_all)

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])

            if i * ncols + j < nattrs:
                attr = attrs[i * ncols + j]
                alpha_value = alpha[i * ncols + j]

                if attr in data.obs:
                    values = data.obs[attr].values
                else:
                    try:
                        pos = data.var_names.get_loc(attr)
                    except KeyError:
                        raise KeyError(f"{attr} is neither in data.obs nor data.var_names!")
                    values = data.X[:, pos].toarray().ravel() if issparse(data.X) else data.X[:, pos]

                if is_numeric_dtype(values):
                    assert apply_to_all or (not restr_obj.contains(attr))
                    values[unsel] = 0 # 0 is good for both integer and float data types

                    img = ax.scatter(
                        x,
                        y,
                        c=values,
                        s=marker_size,
                        marker=".",
                        alpha=alpha_value,
                        edgecolors="none",
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        rasterized=True,
                    )
                    left, bottom, width, height = ax.get_position().bounds
                    rect = [left + width * (1.0 + 0.05), bottom, width * 0.1, height]
                    ax_colorbar = fig.add_axes(rect)
                    fig.colorbar(img, cax=ax_colorbar)

                else:
                    labels = values.astype(str)
                    idx = restr_obj.get_unsatisfied_per_attr(labels, attr) if (not apply_to_all) and restr_obj.contains(attr) else unsel
                    labels[idx] = ""
                    labels = pd.Categorical(labels, categories=natsorted(np.unique(labels)))
                    label_size = labels.categories.size

                    palettes_list = _get_palettes(label_size, with_background=idx.sum() > 0, show_background=show_background) if palettes is None else np.array(palettes.split(","))

                    text_list = []
                    for k, cat in enumerate(labels.categories):
                        idx = labels == cat
                        kwargs = {"marker": ".", "alpha": alpha_value, "edgecolors": "none", "rasterized": True}

                        if legend_loc != "on data":
                            kwargs["label"] = cat
                        else:
                            text_list.append((np.median(x[idx]), np.median(y[idx]), cat))

                        ax.scatter(
                            x[idx],
                            y[idx],
                            c=palettes_list[k],
                            s=marker_size,
                            **kwargs,
                        )

                    if legend_loc == "right margin":
                        legend = ax.legend(
                            loc="center left",
                            bbox_to_anchor=(1, 0.5),
                            frameon=False,
                            fontsize=legend_fontsize,
                            ncol=_get_legend_ncol(label_size, legend_ncol),
                        )
                        for handle in legend.legendHandles:
                            handle.set_sizes([300.0])
                    elif legend_loc == "on data":
                        for px, py, txt in text_list:
                            ax.text(px, py, txt, fontsize=legend_fontsize, horizontalalignment="center", verticalalignment="center")

                ax.set_title(attr)
            else:
                ax.set_frame_on(False)

            if i == nrows - 1:
                ax.set_xlabel(f"{basis}1")
            if j == 0:
                ax.set_ylabel(f"{basis}2")

    return fig if not show else None


def scatter_groups(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attr: str,
    group: str,
    basis: Optional[str] = "umap",
    matkey: Optional[str] = None,
    alpha: Optional[float] = 1.0,
    legend_fontsize: Optional[int] = None,
    show_full: Optional[bool] = True,
    categories: Optional[List[str]] = None,
    palettes: Optional[List[str]] = None,
    legend_ncol: Optional[str] = None,
    cmap: Optional[str] = "YlOrRd",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    subplot_size: Optional[Tuple[float, float]] = (4, 4),
    left: Optional[float] = 0.2,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.4,
    hspace: Optional[float] = 0.15,
    show: Optional[bool] = True,
    **others,
):
    """
    show_full: show the picture with all groups as the first plot.
    categories: if not None, use this to define new groups.
    ### Sample usage:
    ###    fig = plot_scatter_groups(data, 'louvain_labels', 'Individual', 'tsne', nrows = 2, ncols = 4, alpha = 0.5)
    """
    if matkey is not None:
        data.select_matrix(matkey)

    x = data.obsm[f"X_{basis}"][:, 0]
    y = data.obsm[f"X_{basis}"][:, 1]

    basis = _transform_basis(basis)
    marker_size = _get_marker_size(x.size)

    assert group in data.obs
    groups = data.obs[group].values
    if not is_categorical(groups):
        groups = pd.Categorical(groups, categories=natsorted(np.unique(groups)))

    df_g = pd.DataFrame()
    if show_full:
        df_g["All"] = np.ones(data.shape[0], dtype=bool)
    if categories is None:
        for cat in groups.categories:
            df_g[cat] = groups == cat
    else:
        cat_obj = RestrictionParser(categories)
        for key in restr_obj.get_attrs():
            df_g[key] = restr_ojb.get_satisfied_per_attr(groups, key)

    nrows, ncols = _get_nrows_and_ncols(df_g.shape[1], nrows, ncols)
    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=ncols, subplot_size=subplot_size, left=left, bottom=bottom, wspace=wspace, hspace=hspace, squeeze=False)

    if legend_fontsize is None:
        legend_fontsize = rcParams["legend.fontsize"]


    if attr in data.obs:
        values = data.obs[attr].values
    else:
        try:
            pos = data.var_names.get_loc(attr)
        except KeyError:
            raise KeyError(f"{attr} is neither in data.obs nor data.var_names!")
        values = data.X[:, pos].toarray().ravel() if issparse(data.X) else data.X[:, pos]

    is_cat = is_categorical(values)
    if is_cat:
        labels = values
        label_size = labels.categories.size
        palettes = _get_palettes(label_size) if palettes is None else np.array(palettes.split(","))
        legend_ncol = _get_legend_ncol(label_size, legend_ncol)

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])

            gid = i * ncols + j
            if gid < df_g.shape[1]:
                if is_cat:
                    for k, cat in enumerate(labels.categories):
                        idx = np.logical_and(df_g.iloc[:, gid].values, labels == cat)
                        ax.scatter(
                            x[idx],
                            y[idx],
                            c=palettes[k],
                            s=marker_size,
                            marker=".",
                            alpha=alpha,
                            edgecolors="none",
                            label=str(cat),
                            rasterized=True,
                        )

                    legend = ax.legend(
                        loc="center left",
                        bbox_to_anchor=(1, 0.5),
                        frameon=False,
                        fontsize=legend_fontsize,
                        ncol=legend_ncol,
                    )
                    for handle in legend.legendHandles:
                        handle.set_sizes([300.0])
                else:
                    idx_g = df_g.iloc[:, gid].values
                    img = ax.scatter(
                        x[idx_g],
                        y[idx_g],
                        s=marker_size,
                        c=values[idx_g],
                        marker=".",
                        alpha=alpha,
                        edgecolors="none",
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        rasterized=True,
                    )
                    left, bottom, width, height = ax.get_position().bounds
                    rect = [left + width * (1.0 + 0.05), bottom, width * 0.1, height]
                    ax_colorbar = fig.add_axes(rect)
                    fig.colorbar(img, cax=ax_colorbar)

                ax.set_title(str(df_g.columns[gid]))
            else:
                ax.set_frame_on(False)

            if i == nrows - 1:
                ax.set_xlabel(basis + "1")
            if j == 0:
                ax.set_ylabel(basis + "2")

    return fig if not show else None


# def plot_composition(
#     data,
#     cluster,
#     attr,
#     style="frequency",
#     stacked=True,
#     logy=False,
#     subplot_size=(6, 4),
#     left=0.15,
#     bottom=None,
#     wspace=0.3,
#     hspace=None,
#     restrictions=[],
#     **others,
# ):
#     """Generate a composition plot, which shows the percentage of cells from each condition for every cluster.

#     This function is used to generate composition plots, which are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects.

#     Parameters
#     ----------

#     data : `AnnData` or `UnimodalData` or `MultimodalData` object
#         Single cell expression data as an anndata object.
#     cluster : `str`
#         A string represents cluster labels, e.g. louvain_labels.
#     attr: `str`
#         A sample attribute representing the condition, e.g. Donor.
#     style: `str`, optional (default: `frequency`)
#         Composition plot style. Can be either `frequency`, `count`, or 'normalized'. Within each cluster, the `frequency` style show the ratio of cells from each condition over all cells in the cluster, the `count` style just shows the number of cells from each condition, the `normalized` style shows the percentage of cells from the condition in this cluster over the total number of cells from the condition for each condition.
#     stacked: `bool`, optional (default: `True`)
#         If stack the bars from each condition.
#     logy: `bool`, optional (default: `False`)
#         If show the y-axis in log10 scale
#     subplot_size: `tuple`, optional (default: `(6, 4)`)
#         The plot size (width, height) in inches.
#     left: `float`, optional (default: `0.15`)
#         This parameter sets the figure's left margin as a fraction of subplot's width (left * subplot_size[0]).
#     bottom: `float`, optional (default: `0.15`)
#         This parameter sets the figure's bottom margin as a fraction of subplot's height (bottom * subplot_size[1]),
#     wspace: `float`, optional (default: `0.2`)
#         This parameter sets the width between subplots and also the figure's right margin as a fraction of subplot's width (wspace * subplot_size[0]).
#     hspace: `float`, optional (defualt: `0.15`)
#         This parameter sets the height between subplots and also the figure's top margin as a fraction of subplot's height (hspace * subplot_size[1]).
#     restrictions: `list[str]`, optional (default: `[]`)
#         This parameter is used to select a subset of data to plot.

#     Returns
#     -------

#     `Figure` object
#         A `matplotlib.figure.Figure` object containing the composition plot.

#     Examples
#     --------
#     >>> fig = plotting.plot_composition(data, 'louvain_labels', 'Donor', style = 'normalized', stacked = False)
#     """

#     kwargs = set_up_kwargs(subplot_size, left, bottom, wspace, hspace)
#     fig, ax = get_subplot_layouts(**kwargs)

#     restr_obj = RestrictionParser(restrictions)
#     selected = restr_obj.get_satisfied(data)

#     df = pd.crosstab(data.obs.loc[selected, cluster], data.obs.loc[selected, attr])
#     df = df.reindex(
#         index=natsorted(df.index.values), columns=natsorted(df.columns.values)
#     )
#     # df = df.reindex(index = natsorted(df.index.values), columns = ['singlet', 'doublet', 'unknown'])

#     if style == "frequency":
#         df = df.div(df.sum(axis=1), axis=0) * 100.0
#     elif style == "normalized":
#         df = df.div(df.sum(axis=0), axis=1) * 100.0

#     palettes = get_palettes(df.shape[1])

#     rot = None
#     if len(max(df.index.astype(str), key=len)) < 5:
#         rot = 0

#     if logy and not stacked:
#         df_sum = df.sum(axis=1)
#         df_new = df.cumsum(axis=1)
#         df_new = 10 ** df_new.div(df_sum, axis=0).mul(np.log10(df_sum), axis=0)
#         df = df_new.diff(axis=1).fillna(value=df_new.iloc[:, 0:1], axis=1)
#         df.plot(
#             kind="bar",
#             stacked=False,
#             legend=False,
#             logy=True,
#             ylim=(1.01, df_sum.max() * 1.7),
#             color=palettes,
#             rot=rot,
#             ax=ax,
#         )
#     else:
#         df.plot(
#             kind="bar",
#             stacked=stacked,
#             legend=False,
#             logy=logy,
#             color=palettes,
#             rot=rot,
#             ax=ax,
#         )

#     ax.grid(False)
#     ax.set_xlabel("Cluster ID")
#     ax.set_ylabel("Percentage" if style != "count" else "Count")
#     ax.set_title(
#         "AMI = {0:.4f}".format(
#             adjusted_mutual_info_score(
#                 data.obs.loc[selected, cluster],
#                 data.obs.loc[selected, attr],
#                 average_method="arithmetic",
#             )
#         )
#     )
#     ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

#     return fig








# ### Sample usage:
# ###     cg = plot_heatmap(data, 'louvain_labels', ['CD8A', 'CD4', 'CD3G', 'MS4A1', 'NCAM1', 'CD14', 'ITGAX', 'IL3RA', 'CD38', 'CD34', 'PPBP'], use_raw = True, title="markers")
# ###     cg.savefig("heatmap.png", bbox_inches='tight', dpi=600)
# def plot_heatmap(
#     data, cluster, genes, use_raw=False, showzscore=False, title="", **kwargs
# ):
#     sns.set(font_scale=0.35)

#     adata = data.raw if use_raw else data
#     df = pd.DataFrame(adata[:, genes].X.toarray(), index=data.obs.index, columns=genes)
#     if showzscore:
#         df = df.apply(zscore, axis=0)

#     cluster_ids = as_category(data.obs[cluster])
#     idx = cluster_ids.argsort()
#     df = df.iloc[idx, :]  # organize df by category order
#     row_colors = np.zeros(df.shape[0], dtype=object)
#     palettes = get_palettes(cluster_ids.categories.size)

#     cluster_ids = cluster_ids[idx]
#     for k, cat in enumerate(cluster_ids.categories):
#         row_colors[np.isin(cluster_ids, cat)] = palettes[k]

#     cg = sns.clustermap(
#         data=df,
#         row_colors=row_colors,
#         row_cluster=False,
#         col_cluster=True,
#         linewidths=0,
#         yticklabels=[],
#         xticklabels=genes
#     )
#     cg.ax_heatmap.set_ylabel("")
#     # move the colorbar
#     cg.ax_row_dendrogram.set_visible(False)
#     dendro_box = cg.ax_row_dendrogram.get_position()
#     dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
#     dendro_box.x1 = dendro_box.x0 + 0.02
#     cg.cax.set_position(dendro_box)
#     cg.cax.yaxis.set_ticks_position("left")
#     cg.cax.tick_params(labelsize=10)
#     # draw a legend for the cluster groups
#     cg.ax_col_dendrogram.clear()
#     for k, cat in enumerate(cluster_ids.categories):
#         cg.ax_col_dendrogram.bar(0, 0, color=palettes[k], label=cat, linewidth=0)
#     cg.ax_col_dendrogram.legend(loc="center", ncol=15, fontsize=10)
#     cg.ax_col_dendrogram.grid(False)
#     cg.ax_col_dendrogram.set_xticks([])
#     cg.ax_col_dendrogram.set_yticks([])

#     return cg


# ### Sample usage:
# ###     cg = plot_violin_genes(data, 'louvain_labels', ['CD8A', 'CD4', 'CD3G', 'MS4A1', 'NCAM1', 'CD14', 'ITGAX', 'IL3RA', 'CD38', 'CD34', 'PPBP'], use_raw = True, title="markers")
# ###     cg.savefig("heatmap.png", bbox_inches='tight', dpi=600)
# def plot_violin_genes(data, cluster, genes, subplot_size, ylab):
#     nrows, ncols = get_nrows_and_ncols(len(genes), None, None)
#     fig, axes = get_subplot_layouts(
#         nrows=nrows, ncols=ncols, subplot_size=subplot_size, hspace=0.3, wspace=0.1, squeeze=False
#     )
#     expr_mat = data[:, genes].X.toarray()
#     df = pd.DataFrame(data=expr_mat, columns=genes)
#     df.insert(0, "label", data.obs[cluster].values)
#     for i in range(nrows):
#         for j in range(ncols):
#             ax = axes[i, j]
#             idx = i * ncols + j
#             if idx < len(genes):
#                 sns.violinplot(
#                     x="label", y=genes[idx], data=df, inner=None, linewidth=0, ax=ax, cut=0
#                 )
#                 sns.stripplot(x="label", y=genes[idx], data=df, ax=ax, size=2, color="k")
#                 ax.set_xlabel("")
#                 ax.set_ylabel("")
#                 ax.set_title(genes[idx])
#             else:
#                 ax.set_frame_on(False)
#     plt.figtext(0.02, 0.5, ylab, rotation="vertical", fontsize="xx-large")
#     return fig
