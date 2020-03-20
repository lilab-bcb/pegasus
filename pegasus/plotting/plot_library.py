import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import zscore
from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted
from matplotlib import rcParams
from pandas.api.types import is_numeric_dtype
from collections import namedtuple

from typing import List

from .plot_utils import get_palettes, transform_basis


Restriction = namedtuple("negation", "values")


class RestrictionParser:
    def __init__(self, restrictions: List[str]):
        self.restrs = {}
        for restr_str in restrictions:
            attr, value_str = restr_str.split(":")
            negation = False
            if value_str[0] == "~":
                negation = True
                value_str = value_str[1:]
            self.restr[attr] = Restriction(
                negation=negation, values=value_str.split(",")
            )

    def contains(self, attr: str) -> bool:
        return attr in self.restrs

    def get_attrs(self) -> List[str]:
        return self.restrs.keys()

    def get_satisfied(self, data: "AnnData") -> List[bool]:
        selected = np.ones(data.shape[0], dtype=bool)
        for attr, restr in self.restrs.items():
            labels = data.obs[attr].astype(str)
            if restr.negation:
                selected = selected & (~np.isin(labels, restr.values))
            else:
                selected = selected & np.isin(labels, restr.values)
        return selected

    def get_unsatisfied(self, data: "AnnData", apply_to_all: bool) -> List[bool]:
        unsel = np.zeros(data.shape[0], dtype=bool)
        if apply_to_all:
            for attr, restr in self.restrs.items():
                labels = data.obs[attr].astype(str)
                if restr.negation:
                    unsel = unsel | np.isin(labels, restr.values)
                else:
                    unsel = unsel | (~np.isin(labels, restr.values))
        return unsel

    def get_satisfied_per_attr(self, labels: List[str], attr: str) -> List[bool]:
        one_restr = self.restrs[attr]
        if one_restr.negation:
            return ~np.isin(labels, rest_vec)
        else:
            return np.isin(labels, rest_vec)

    def get_unsatisfied_per_attr(self, labels: List[str], attr: str) -> List[bool]:
        one_restr = self.restrs[attr]
        if one_restr.negation:
            return np.isin(labels, rest_vec)
        else:
            return ~np.isin(labels, rest_vec)


def get_nrows_and_ncols(num_figs, nrows, ncols):
    if nrows is None and ncols is None:
        nrows = int(np.sqrt(num_figs))
        ncols = (num_figs // nrows) + (num_figs % nrows > 0)
    elif nrows is None:
        nrows = (num_figs // ncols) + (num_figs % ncols > 0)
    elif ncols is None:
        ncols = (num_figs // nrows) + (num_figs % nrows > 0)

    return nrows, ncols


def set_up_kwargs(subplot_size, left, bottom, wspace, hspace):
    kwargs = {}
    if subplot_size is not None:
        kwargs["subplot_size"] = subplot_size
    if left is not None:
        kwargs["left"] = left
    if bottom is not None:
        kwargs["bottom"] = bottom
    if wspace is not None:
        kwargs["wspace"] = wspace
    if hspace is not None:
        kwargs["hspace"] = hspace

    return kwargs


def get_subplot_layouts(
    nrows=1,
    ncols=1,
    subplot_size=(4, 4),
    left=0.2,
    bottom=0.15,
    wspace=0.4,
    hspace=0.15,
    squeeze=True,
    sharex=True,
    sharey=True,
    frameon=False,
):
    left_margin = left * subplot_size[0]
    bottom_margin = bottom * subplot_size[1]
    right_space = wspace * subplot_size[0]
    top_space = hspace * subplot_size[1]

    figsize = (
        left_margin + subplot_size[0] * (1.0 + wspace) * ncols,
        bottom_margin + subplot_size[1] * (1.0 + hspace) * nrows,
    )
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        squeeze=squeeze,
        sharex=sharex,
        sharey=sharey,
        frameon=frameon,
    )

    fig.subplots_adjust(
        left=left_margin / figsize[0],
        bottom=bottom_margin / figsize[1],
        right=1.0 - right_space / figsize[0],
        top=1.0 - top_space / figsize[1],
        wspace=wspace,
        hspace=hspace,
    )

    return fig, axes


def get_legend_ncol(label_size):
    return 1 if label_size <= 14 else (2 if label_size <= 30 else 3)


def as_category(labels):
    if isinstance(labels, pd.Series):
        labels = labels.values
    # if 'singlet' in np.unique(labels):
    #     return pd.Categorical(labels, categories = ['', 'singlet', 'doublet', 'unknown'])
    return pd.Categorical(labels, categories=natsorted(np.unique(labels)))


def get_marker_size(nsamples):
    return (240000.0 if nsamples > 300000 else 120000.0) / nsamples


def plot_composition(
    data,
    cluster,
    attr,
    style="frequency",
    stacked=True,
    logy=False,
    subplot_size=(6, 4),
    left=0.15,
    bottom=None,
    wspace=0.3,
    hspace=None,
    restrictions=[],
    **others,
):
    """Generate a composition plot, which shows the percentage of cells from each condition for every cluster.
    
    This function is used to generate composition plots, which are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects.

    Parameters
    ----------

    data : `anndata` object
        Single cell expression data as an anndata object.
    cluster : `str`
        A string represents cluster labels, e.g. louvain_labels.
    attr: `str`
        A sample attribute representing the condition, e.g. Donor.
    style: `str`, optional (default: `frequency`)
        Composition plot style. Can be either `frequency`, `count`, or 'normalized'. Within each cluster, the `frequency` style show the ratio of cells from each condition over all cells in the cluster, the `count` style just shows the number of cells from each condition, the `normalized` style shows the percentage of cells from the condition in this cluster over the total number of cells from the condition for each condition. 
    stacked: `bool`, optional (default: `True`)
        If stack the bars from each condition.
    logy: `bool`, optional (default: `False`)
        If show the y-axis in log10 scale
    subplot_size: `tuple`, optional (default: `(6, 4)`)
        The plot size (width, height) in inches.
    left: `float`, optional (default: `0.15`)
        This parameter sets the figure's left margin as a fraction of subplot's width (left * subplot_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of subplot's height (bottom * subplot_size[1]),
    wspace: `float`, optional (default: `0.2`)
        This parameter sets the width between subplots and also the figure's right margin as a fraction of subplot's width (wspace * subplot_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between subplots and also the figure's top margin as a fraction of subplot's height (hspace * subplot_size[1]).
    restrictions: `list[str]`, optional (default: `[]`)
        This parameter is used to select a subset of data to plot.
    
    Returns
    -------

    `Figure` object
        A `matplotlib.figure.Figure` object containing the composition plot.

    Examples
    --------
    >>> fig = plotting.plot_composition(data, 'louvain_labels', 'Donor', style = 'normalized', stacked = False)
    """

    kwargs = set_up_kwargs(subplot_size, left, bottom, wspace, hspace)
    fig, ax = get_subplot_layouts(**kwargs)

    restr_obj = RestrictionParser(restrictions)
    selected = restr_obj.get_satisfied(data)

    df = pd.crosstab(data.obs.loc[selected, cluster], data.obs.loc[selected, attr])
    df = df.reindex(
        index=natsorted(df.index.values), columns=natsorted(df.columns.values)
    )
    # df = df.reindex(index = natsorted(df.index.values), columns = ['singlet', 'doublet', 'unknown'])

    if style == "frequency":
        df = df.div(df.sum(axis=1), axis=0) * 100.0
    elif style == "normalized":
        df = df.div(df.sum(axis=0), axis=1) * 100.0

    palettes = get_palettes(df.shape[1])

    rot = None
    if len(max(df.index.astype(str), key=len)) < 5:
        rot = 0

    if logy and not stacked:
        df_sum = df.sum(axis=1)
        df_new = df.cumsum(axis=1)
        df_new = 10 ** df_new.div(df_sum, axis=0).mul(np.log10(df_sum), axis=0)
        df = df_new.diff(axis=1).fillna(value=df_new.iloc[:, 0:1], axis=1)
        df.plot(
            kind="bar",
            stacked=False,
            legend=False,
            logy=True,
            ylim=(1.01, df_sum.max() * 1.7),
            color=palettes,
            rot=rot,
            ax=ax,
        )
    else:
        df.plot(
            kind="bar",
            stacked=stacked,
            legend=False,
            logy=logy,
            color=palettes,
            rot=rot,
            ax=ax,
        )

    ax.grid(False)
    ax.set_xlabel("Cluster ID")
    ax.set_ylabel("Percentage" if style != "count" else "Count")
    ax.set_title(
        "AMI = {0:.4f}".format(
            adjusted_mutual_info_score(
                data.obs.loc[selected, cluster],
                data.obs.loc[selected, attr],
                average_method="arithmetic",
            )
        )
    )
    ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

    return fig


### Sample usage:
###    fig = plot_scatter(data, 'tsne', ['louvain_labels', 'hdbscan_labels_soft'], nrows = 1, ncols = 2, alpha = 0.5)
def plot_scatter(
    data,
    basis,
    attrs,
    restrictions=[],
    nrows=None,
    ncols=None,
    subplot_size=(4, 4),
    left=None,
    bottom=None,
    wspace=None,
    hspace=None,
    alpha=None,
    legend_fontsize=None,
    apply_to_all=True,
    show_background=False,
    vmin=0.0,
    vmax=1.0,
    **others,
):
    df = pd.DataFrame(
        data.obsm["X_" + basis][:, 0:2], columns=[basis + c for c in ["1", "2"]]
    )
    basis = transform_basis(basis)

    nattrs = len(attrs)
    nrows, ncols = get_nrows_and_ncols(nattrs, nrows, ncols)
    marker_size = get_marker_size(df.shape[0])

    kwargs = set_up_kwargs(subplot_size, left, bottom, wspace, hspace)
    fig, axes = get_subplot_layouts(nrows=nrows, ncols=ncols, squeeze=False, **kwargs)

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
                alpha_value = alpha[i * ncols + j] if isinstance(alpha, list) else alpha

                if is_numeric_dtype(data.obs[attr]):
                    values = data.obs[attr].values
                    assert apply_to_all or (not restr_obj.contains(attr))
                    values[unsel] = 0.0

                    img = ax.scatter(
                        df.iloc[:, 0],
                        df.iloc[:, 1],
                        c=values,
                        s=marker_size,
                        marker=".",
                        alpha=alpha_value,
                        edgecolors="none",
                        cmap="viridis",
                        rasterized=True,
                    )
                    img.set_clim(vmin, vmax)
                    left, bottom, width, height = ax.get_position().bounds
                    rect = [left + width * (1.0 + 0.05), bottom, width * 0.1, height]
                    ax_colorbar = fig.add_axes(rect)
                    fig.colorbar(img, cax=ax_colorbar)

                else:
                    labels = data.obs[attr].astype(str)
                    if (not apply_to_all) and restr_obj.contains(attr):
                        idx = restr_obj.get_unsatisfied_per_attr(labels, attr)
                    else:
                        idx = unsel
                    labels[idx] = ""
                    labels = as_category(labels)
                    label_size = labels.categories.size

                    if others["palettes"] is not None:
                        palettes = np.array(others["palettes"].split(","))
                    else:
                        palettes = get_palettes(
                            label_size,
                            with_background=idx.sum() > 0,
                            show_background=show_background,
                        )

                    for k, cat in enumerate(labels.categories):
                        idx = np.isin(labels, cat)
                        ax.scatter(
                            df.iloc[idx, 0],
                            df.iloc[idx, 1],
                            c=palettes[k],
                            s=marker_size,
                            marker=".",
                            alpha=alpha_value,
                            edgecolors="none",
                            label=cat,
                            rasterized=True,
                        )

                    legend = ax.legend(
                        loc="center left",
                        bbox_to_anchor=(1, 0.5),
                        frameon=False,
                        fontsize=legend_fontsize,
                        ncol=get_legend_ncol(label_size),
                    )
                    for handle in legend.legendHandles:
                        handle.set_sizes([300.0])

                ax.set_title(attr)
            else:
                ax.set_frame_on(False)

            if i == nrows - 1:
                ax.set_xlabel(basis + "1")
            if j == 0:
                ax.set_ylabel(basis + "2")

    return fig


### Sample usage:
###    fig = plot_scatter_groups(data, 'tsne', 'louvain_labels', 'Individual', nrows = 2, ncols = 4, alpha = 0.5)
def plot_scatter_groups(
    data,
    basis,
    cluster,
    group,
    restrictions=[],
    nrows=None,
    ncols=None,
    subplot_size=(4, 4),
    left=None,
    bottom=None,
    wspace=None,
    hspace=None,
    alpha=None,
    legend_fontsize=None,
    showall=True,
    **others,
):
    df = pd.DataFrame(
        data.obsm["X_" + basis][:, 0:2], columns=[basis + c for c in ["1", "2"]]
    )
    basis = transform_basis(basis)

    marker_size = get_marker_size(df.shape[0])

    labels = as_category(data.obs[cluster])
    label_size = labels.categories.size
    palettes = get_palettes(label_size)
    legend_ncol = get_legend_ncol(label_size)

    assert group in data.obs
    groups = as_category(data.obs[group])
    df_g = pd.DataFrame()

    if showall:
        df_g["All"] = np.ones(data.shape[0], dtype=bool)
    if len(restrictions) == 0:
        for cat in groups.categories:
            df_g[cat] = np.isin(groups, cat)
    else:
        restr_obj = RestrictionParser(restrictions)
        for key in restr_obj.get_attrs():
            df_g[key] = restr_ojb.get_satisfied_per_attr(groups, key)

    nrows, ncols = get_nrows_and_ncols(df_g.shape[1], nrows, ncols)

    kwargs = set_up_kwargs(subplot_size, left, bottom, wspace, hspace)
    fig, axes = get_subplot_layouts(nrows=nrows, ncols=ncols, squeeze=False, **kwargs)

    if legend_fontsize is None:
        legend_fontsize = rcParams["legend.fontsize"]

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])

            gid = i * ncols + j
            if gid < df_g.shape[1]:
                for k, cat in enumerate(labels.categories):
                    idx = np.logical_and(df_g.iloc[:, gid].values, np.isin(labels, cat))
                    ax.scatter(
                        df.iloc[idx, 0],
                        df.iloc[idx, 1],
                        c=palettes[k],
                        s=marker_size,
                        marker=".",
                        alpha=alpha,
                        edgecolors="none",
                        label=cat,
                        rasterized=True,
                    )

                ax.set_title(df_g.columns[gid])
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
                ax.set_frame_on(False)

            if i == nrows - 1:
                ax.set_xlabel(basis + "1")
            if j == 0:
                ax.set_ylabel(basis + "2")

    return fig


### Sample usage:
###    fig = plot_scatter_genes(data, 'tsne', ['CD8A', 'CD4', 'CD3G', 'MS4A1', 'NCAM1', 'CD14', 'ITGAX', 'IL3RA', 'CD38', 'CD34', 'PPBP'])
def plot_scatter_genes(
    data,
    basis,
    genes,
    nrows=None,
    ncols=None,
    subplot_size=(4, 4),
    left=None,
    bottom=None,
    wspace=0.3,
    hspace=None,
    alpha=None,
    use_raw=False,
    **others,
):
    df = pd.DataFrame(
        data.obsm["X_" + basis][:, 0:2], columns=[basis + c for c in ["1", "2"]]
    )
    basis = transform_basis(basis)

    ngenes = len(genes)
    nrows, ncols = get_nrows_and_ncols(ngenes, nrows, ncols)
    marker_size = 240000.0 / df.shape[0]

    kwargs = set_up_kwargs(subplot_size, left, bottom, wspace, hspace)
    fig, axes = get_subplot_layouts(nrows=nrows, ncols=ncols, squeeze=False, **kwargs)

    X = data.raw[:, genes].X.toarray() if use_raw else data[:, genes].X.toarray()

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])

            if i * ncols + j < ngenes:
                gene_id = i * ncols + j
                img = ax.scatter(
                    df.iloc[:, 0],
                    df.iloc[:, 1],
                    s=marker_size,
                    c=X[:, gene_id],
                    marker=".",
                    alpha=alpha,
                    edgecolors="none",
                    cmap="viridis",
                    rasterized=True,
                )

                left, bottom, width, height = ax.get_position().bounds
                rect = [left + width * (1.0 + 0.05), bottom, width * 0.1, height]
                ax_colorbar = fig.add_axes(rect)
                fig.colorbar(img, cax=ax_colorbar)

                ax.set_title(genes[gene_id])
            else:
                ax.set_frame_on(False)

            if i == nrows - 1:
                ax.set_xlabel(basis + "1")
            if j == 0:
                ax.set_ylabel(basis + "2")

    return fig


### Sample usage:
###    fig = plot_scatter_gene_groups(data, 'tsne', 'CD8A', 'Individual', nrows = 3, ncols = 3)
def plot_scatter_gene_groups(
    data,
    basis,
    gene,
    group,
    nrows=None,
    ncols=None,
    subplot_size=(4, 4),
    left=None,
    bottom=None,
    wspace=0.3,
    hspace=None,
    alpha=None,
    use_raw=False,
    **others,
):
    df = pd.DataFrame(
        data.obsm["X_" + basis][:, 0:2], columns=[basis + c for c in ["1", "2"]]
    )
    basis = transform_basis(basis)

    marker_size = 240000.0 / df.shape[0]
    groups = as_category(data.obs[group])
    ngroup = groups.categories.size
    nrows, ncols = get_nrows_and_ncols(ngroup + 1, nrows, ncols)

    kwargs = set_up_kwargs(subplot_size, left, bottom, wspace, hspace)
    fig, axes = get_subplot_layouts(nrows=nrows, ncols=ncols, squeeze=False, **kwargs)

    X = (data.raw[:, gene].X if use_raw else data[:, gene].X).toarray()[:, 0]
    vmax = X.max()

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])

            gid = i * ncols + j
            if gid <= ngroup:
                if gid == 0:
                    idx_g = np.ones(groups.shape[0], dtype=bool)
                else:
                    idx_g = np.isin(groups, groups.categories[gid - 1])

                img = ax.scatter(
                    df.iloc[idx_g, 0],
                    df.iloc[idx_g, 1],
                    s=marker_size,
                    c=X[idx_g],
                    marker=".",
                    alpha=alpha,
                    edgecolors="none",
                    cmap="viridis",
                    vmin=0.0,
                    vmax=vmax,
                    rasterized=True,
                )

                left, bottom, width, height = ax.get_position().bounds
                rect = [left + width * (1.0 + 0.05), bottom, width * 0.1, height]
                ax_colorbar = fig.add_axes(rect)
                fig.colorbar(img, cax=ax_colorbar)

                ax.set_title("All" if gid == 0 else str(groups.categories[gid - 1]))
            else:
                ax.set_frame_on(False)

            if i == nrows - 1:
                ax.set_xlabel(basis + "1")
            if j == 0:
                ax.set_ylabel(basis + "2")

    return fig


### Sample usage:
###     cg = plot_heatmap(data, 'louvain_labels', ['CD8A', 'CD4', 'CD3G', 'MS4A1', 'NCAM1', 'CD14', 'ITGAX', 'IL3RA', 'CD38', 'CD34', 'PPBP'], use_raw = True, title="markers")
###     cg.savefig("heatmap.png", bbox_inches='tight', dpi=600)
def plot_heatmap(
    data, cluster, genes, use_raw=False, showzscore=False, title="", **kwargs
):
    sns.set(font_scale=0.35)

    adata = data.raw if use_raw else data
    df = pd.DataFrame(adata[:, genes].X.toarray(), index=data.obs.index, columns=genes)
    if showzscore:
        df = df.apply(zscore, axis=0)

    cluster_ids = as_category(data.obs[cluster])
    idx = cluster_ids.argsort()
    df = df.iloc[idx, :]  # organize df by category order
    row_colors = np.zeros(df.shape[0], dtype=object)
    palettes = get_palettes(cluster_ids.categories.size)

    cluster_ids = cluster_ids[idx]
    for k, cat in enumerate(cluster_ids.categories):
        row_colors[np.isin(cluster_ids, cat)] = palettes[k]

    cg = sns.clustermap(
        data=df,
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=True,
        linewidths=0,
        yticklabels=[],
        xticklabels=genes
    )
    cg.ax_heatmap.set_ylabel("")
    # move the colorbar
    cg.ax_row_dendrogram.set_visible(False)
    dendro_box = cg.ax_row_dendrogram.get_position()
    dendro_box.x0 = (dendro_box.x0 + 2 * dendro_box.x1) / 3
    dendro_box.x1 = dendro_box.x0 + 0.02
    cg.cax.set_position(dendro_box)
    cg.cax.yaxis.set_ticks_position("left")
    cg.cax.tick_params(labelsize=10)
    # draw a legend for the cluster groups
    cg.ax_col_dendrogram.clear()
    for k, cat in enumerate(cluster_ids.categories):
        cg.ax_col_dendrogram.bar(0, 0, color=palettes[k], label=cat, linewidth=0)
    cg.ax_col_dendrogram.legend(loc="center", ncol=15, fontsize=10)
    cg.ax_col_dendrogram.grid(False)
    cg.ax_col_dendrogram.set_xticks([])
    cg.ax_col_dendrogram.set_yticks([])

    return cg


### Sample usage:
###     cg = plot_violin_genes(data, 'louvain_labels', ['CD8A', 'CD4', 'CD3G', 'MS4A1', 'NCAM1', 'CD14', 'ITGAX', 'IL3RA', 'CD38', 'CD34', 'PPBP'], use_raw = True, title="markers")
###     cg.savefig("heatmap.png", bbox_inches='tight', dpi=600)
def plot_violin_genes(data, cluster, genes, subplot_size, ylab):
    nrows, ncols = get_nrows_and_ncols(len(genes), None, None)
    fig, axes = get_subplot_layouts(
        nrows=nrows, ncols=ncols, subplot_size=subplot_size, hspace=0.3, wspace=0.1, squeeze=False
    )
    expr_mat = data[:, genes].X.toarray()
    df = pd.DataFrame(data=expr_mat, columns=genes)
    df.insert(0, "label", data.obs[cluster].values)
    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            idx = i * ncols + j
            if idx < len(genes):
                sns.violinplot(
                    x="label", y=genes[idx], data=df, inner=None, linewidth=0, ax=ax, cut=0
                )
                sns.stripplot(x="label", y=genes[idx], data=df, ax=ax, size=2, color="k")
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.set_title(genes[idx])
            else:
                ax.set_frame_on(False)
    plt.figtext(0.02, 0.5, ylab, rotation="vertical", fontsize="xx-large")
    return fig
