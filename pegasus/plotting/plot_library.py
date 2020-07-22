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

from typing import List, Tuple, Union, Optional, Callable

import logging
logger = logging.getLogger(__name__)

from pegasus.tools import X_from_rep
from .plot_utils import _transform_basis, _get_nrows_and_ncols, _get_marker_size, _get_dot_size, _get_subplot_layouts, _get_legend_ncol, _get_palettes, RestrictionParser


def scatter(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attrs: Union[str, List[str]],
    basis: Optional[str] = "umap",
    matkey: Optional[str] = None,
    alpha: Optional[Union[float, List[float]]] = 1.0,
    legend_loc: Optional[str] = "right margin",
    legend_fontsize: Optional[int] = None,
    restrictions: Optional[List[str]] = None,
    apply_to_all: Optional[bool] = True,
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
    """Generate scatter plots for different attributes

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
       Use current selected modality in data.
    attrs: ``str`` or ``List[str]``
        Color scatter plots by attrs. Each attribute in attrs should be one key in data.obs or data.var_names (e.g. one gene). If one attribute is categorical, a palette will be used to color each category separately. Otherwise, a color map will be used.
    basis: ``str``, optional, default: ``umap``
        Basis to be used to generate scatter plots. Can be either 'umap', 'tsne', 'fitsne', 'fle', 'net_tsne', 'net_fitsne', 'net_umap' or 'net_fle'.
    matkey: ``str``, optional, default: None
        If matkey is set, select matrix with matkey as keyword in the current modality. Only works for MultimodalData or UnimodalData objects.
    alpha: ``float`` or ``List[float], optional, default: ``1.0``
        Alpha value for blending, from 0.0 (transparent) to 1.0 (opaque). If this is a list, the length must match attrs, which means we set a separate alpha value for each attribute.
    legend_loc: ``str``, optional, default: ``right margin``
        Legend location. Can be either "right margin" or "on data".
    legend_fontsize: ``int``, optional, default: None
        Legend fontsize. If None, use matplotlib.rcParams["legend.fontsize"].
    restrictions: ``List[str]``, optional, default: None

    Returns
    -------
    ``None``

    Update ``data.obs``:

        * ``n_genes``: Total number of genes for each cell.
        * ``n_counts``: Total number of counts for each cell.
        * ``percent_mito``: Percent of mitochondrial genes for each cell.
        * ``passed_qc``: Boolean type indicating if a cell passes the QC process based on the QC metrics.
        * ``demux_type``: this column might be deleted if select_singlets is on.

    Examples
    --------
    >>> pg.qc_metrics(data)
    """



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

    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()

    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
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

    # Reset current matrix if needed.
    if (not isinstance(data, anndata.AnnData)):
        if cur_matkey != data.current_matrix():
            data.select_matrix(cur_matkey)

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
    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
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

    if not isinstance(data, anndata.AnnData):
        if cur_matkey != data.current_matrix():
            data.select_matrix(cur_matkey)

    return fig if not show else None


def compo_plot(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    xattr: str,
    yattr: str,
    style: Optional[str] = "frequency",
    restrictions: Optional[List[str]] = None,
    xlabel: Optional[str] = None,
    subplot_size: Optional[Tuple[float, float]] = (6, 4),
    left: Optional[float] = 0.15,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.3,
    hspace: Optional[float] = 0.15,
    show: Optional[bool] = True,
    **others,
):
    """Generate a composition plot, which shows the percentage of cells from each condition for every cluster.

    This function is used to generate composition plots, which are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects.

    Parameters
    ----------

    data : `AnnData` or `UnimodalData` or `MultimodalData` object
        Single cell expression data.
    xattr : `str`
        A string for the attribute used in x-axis, e.g. Donor.
    yattr: `str`
        A string for the attribute used in y-axis, e.g. Cell type.
    style: `str`, optional (default: `frequency`)
        Composition plot style. Can be either `frequency`, or 'normalized'. Within each cluster, the `frequency` style show the percentage of cells from each yattr over all cells in the xattr (stacked), the `normalized` style shows the percentage of cells from yattr in the xattr over all of cells from yattr for each yattr (not stacked).
    restrictions: `list[str]`, optional (default: None)
        This parameter is used to select a subset of data to plot.
    xlabel: `str`, optional (default None)
        X-axis label. If None, use xattr
    subplot_size: `tuple`, optional (default: `(6, 4)`)
        The plot size (width, height) in inches.
    left: `float`, optional (default: `0.15`)
        This parameter sets the figure's left margin as a fraction of subplot's width (left * subplot_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of subplot's height (bottom * subplot_size[1]).
    wspace: `float`, optional (default: `0.3`)
        This parameter sets the width between subplots and also the figure's right margin as a fraction of subplot's width (wspace * subplot_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between subplots and also the figure's top margin as a fraction of subplot's height (hspace * subplot_size[1]).

    Returns
    -------

    `Figure` object
        A `matplotlib.figure.Figure` object containing the composition plot if show == False

    Examples
    --------
    >>> fig = pg.compo_plot(data, 'louvain_labels', 'Donor', style = 'normalized')
    """
    if xlabel is None:
        xlabel = xattr

    fig, ax = _get_subplot_layouts(nrows=1, ncols=1, subplot_size=subplot_size, left=left, bottom=bottom, wspace=wspace, hspace=hspace)

    restr_obj = RestrictionParser(restrictions)
    selected = restr_obj.get_satisfied(data)

    df = pd.crosstab(data.obs.loc[selected, xattr], data.obs.loc[selected, yattr])
    df = df.reindex(
        index=natsorted(df.index.values), columns=natsorted(df.columns.values)
    )

    if style == "frequency":
        df = df.div(df.sum(axis=1), axis=0) * 100.0
    else:
        assert style == "normalized"
        df = df.div(df.sum(axis=0), axis=1) * 100.0

    palettes = _get_palettes(df.shape[1])

    rot = None
    if len(max(df.index.astype(str), key=len)) < 5:
        rot = 0

    df.plot(
        kind = "bar",
        stacked = style == "frequency",
        legend = False,
        color = palettes,
        rot = rot,
        ax = ax,
    )

    ax.grid(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Percentage")
    ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

    return fig if not show else None


def violin(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    keys: Union[str, List[str]],
    groupby: str,
    matkey: Optional[str] = None,
    stripplot: bool = False,
    scale: str = 'width',
    jitter: Union[float, bool] = False,
    subplot_size: Optional[Tuple[float, float]] = (8, 0.5),
    left: Optional[float] = 0.15,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.1,
    ylabel: Optional[str] = None,
    show: Optional[bool] = True,
    **others,
    ):
    """
    Generate a stacked violin plot.

    Parameters
    ----------
    data: ``AnnData`` or ``MultimodalData`` or ``UnimodalData`` object
        Single-cell expression data.
    keys: ``str`` or ``List[str]``
        Cell attributes or features to plot.
        Cell attributes must exist in ``data.obs`` and must be numeric.
        Features must exist in ``data.var``.
    groupby: ``str``
        Cell attribute to group data points.
    matkey: ``str``, optional, default: ``None``
        If matkey is set, select matrix with matkey as keyword in the current modality. Only works for MultimodalData or UnimodalData objects.
    stripplot: ``bool``, optional, default: ``False``
        Attach a stripplot to the violinplot or not.
    scale: ``str``, optional, default: ``width``
        The method used to scale the width of each violin:
            - If ``width``, each violin will have the same width.
            - If ``area``, each violin will have the same area.
            - If ``count``, the width of the violins will be scaled by the number of observations in that bin.
    jitter: ``float`` or ``bool``, optional, default: ``False``
        Amount of jitter (only along the categorical axis) to apply to stripplot. This is used only when ``stripplot`` is set to ``True``.
        This can be useful when you have many points and they overlap, so that it is easier to see the distribution. You can specify the amount of jitter (half the width of the uniform random variable support), or just use ``True`` for a good default.
    subplot_size: ``Tuple[float, float]``, optional, default: ``(10, 1)``
        The size (width, height) in inches of each violin subplot.
    left: ``float``, optional, default: ``0.15``
        This parameter sets the figure's left margin as a fraction of subplot's width (left * subplot_size[0]).
    bottom: ``float``, optional, default: ``0.15``
        This parameter sets the figure's bottom margin as a fraction of subplot's height (bottom * subplot_size[1]).
    wspace: ``float``, optional, default: ``0.1``
        This parameter sets the width between subplots and also the figure's right margin as a fraction of subplot's width (wspace * subplot_size[0]).
    ylabel: ``str``, optional, default: ``None``
        Y-axis label. No label to show if ``None``.
    show: ``bool``, optional, default: ``True``
        Return a ``Figure`` object if ``False``; return ``None`` otherwise.
    others
        Are passed to ``seaborn.violinplot``.

    Returns
    -------
    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``show == False``

    Examples
    --------
    >>> pg.violin(data, keys=['CD14', 'TRAC', 'CD34'], groupby='louvain_labels')
    """
    if not is_list_like(keys):
        keys = [keys]

    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
        data.select_matrix(matkey)

    nrows, ncols = (len(keys), 1)

    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=ncols, subplot_size=subplot_size, left=left, bottom=bottom, wspace=wspace, hspace=0, squeeze=False, sharey=False)

    obs_keys = []
    genes = []

    for key in keys:
        if key in data.obs:
            assert is_numeric_dtype(data.obs[key])
            obs_keys.append(key)
        else:
            genes.append(key)

    df_list = [pd.DataFrame({"label": data.obs[groupby].values})]
    if len(obs_keys) > 0:
        df_list.append(data.obs[obs_keys])
    if len(genes) > 0:
        expr_mat = data[:, genes].X.toarray()
        df_list.append(pd.DataFrame(data=expr_mat, columns=genes))
    df = pd.concat(df_list, axis = 1)

    for i in range(nrows):
        ax = axes[i, 0]
        if stripplot:
            sns.stripplot(x="label", y=keys[i], data=df, ax=ax, size=1, color="k", jitter=jitter)
        sns.violinplot(x="label", y=keys[i], data=df, inner=None, linewidth=1, ax=ax, cut=0, scale=scale, **others)
        ax.grid(False)
        if i < nrows - 1:
            ax.set_xlabel("")
        else:
            ax.set_xlabel(groupby)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylabel(keys[i], labelpad=8, rotation=0, horizontalalignment='right', fontsize='medium')
        ax.tick_params(axis='y', right=True, left=False, labelright=True, labelleft=False, labelsize='small')

    if ylabel is not None:
        plt.figtext(0.02, 0.5, ylabel, rotation="vertical", fontsize="xx-large")

    # Reset current matrix if needed.
    if not isinstance(data, anndata.AnnData):
        if data.current_matrix() != cur_matkey:
            data.select_matrix(cur_matkey)

    return fig if not show else None


def heatmap(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    cluster: str,
    genes: Union[str, List[str]],
    matkey: Optional[str] = None,
    on_average: bool = True,
    row_cluster: bool = False,
    col_cluster: bool = True,
    show: bool = True,
    **kwargs,
):
    """
    Generate a heatmap.

    Parameters
    -----------

    data: ``AnnData`` or ``MultimodalData`` or ``UnimodalData`` object
        Single-cell expression data.
    cluster: ``str``
        Cell attribute to plot.
    genes: ``str`` or ``List[str]``
        Features to plot.
    matkey: ``str``, optional, default: ``None``
        If matkey is set, select matrix with matkey as keyword in the current modality. Only works for MultimodalData or UnimodalData objects.
    on_average: ``bool``, optional, default: ``True``
        If ``True``, plot cluster average gene expression (i.e. show a Matrixplot); otherwise, plot a general heatmap.
    row_cluster: ``bool``, optional, default: ``False``
    col_cluster: ``bool``, optional, default: ``True``
    cmap: ``str``, optional, default: ``Reds``
        Color map for plotting. See `colormap documentation`_ for a detailed list.
    show: ``bool``, optional, default: ``True``
        Return a ``Figure`` object if ``False``; return ``None`` otherwise.
    kwargs
        Are passed to ``seaborn.heatmap``.

    .. _colormap documentation: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    Returns
    -------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``show == False``

    Examples
    --------
    >>> pg.heatmap(data, cluster='louvain_labels', genes=['CD14', 'TRAC', 'CD34'])

    """
    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
        data.select_matrix(matkey)

    df = pd.DataFrame(data[:, genes].X.toarray(), index=data.obs.index, columns=genes)
    df['cluster_name'] = data.obs[cluster]

    if on_average:
        if not 'cmap' in kwargs.keys():
            kwargs['cmap'] = 'viridis'
        df = df.groupby('cluster_name').mean()
        cluster_ids = df.index
    else:
        cluster_ids = pd.Categorical(data.obs[cluster])
        idx = cluster_ids.argsort()
        df = df.iloc[idx, :]  # organize df by category order
        df.drop(columns=['cluster_name'], inplace=True)

    if not on_average:
        row_colors = np.zeros(df.shape[0], dtype=object)
        palettes = _get_palettes(cluster_ids.categories.size)
        cluster_ids = cluster_ids[idx]
        for k, cat in enumerate(cluster_ids.categories):
            row_colors[np.isin(cluster_ids, cat)] = palettes[k]

    cg = sns.clustermap(
        data=df,
        row_colors=row_colors if not on_average else None,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        linewidths=0,
        yticklabels=cluster_ids if on_average else [],
        xticklabels=genes,
        **kwargs,
    )
    cg.ax_heatmap.set_ylabel("")

    if row_cluster:
        cg.ax_heatmap.yaxis.tick_right()
    else:
        cg.ax_heatmap.yaxis.tick_left()

    cg.ax_row_dendrogram.set_visible(row_cluster)
    cg.cax.tick_params(labelsize=10)
    # move the colorbar if needed
    if not (row_cluster and on_average):
        color_box = cg.ax_heatmap.get_position()
        color_box.x0 = color_box.x1 + 0.04
        color_box.x1 = color_box.x0 + 0.02
        cg.cax.set_position(color_box)
        cg.cax.yaxis.set_ticks_position("right")
    # draw a legend for the cluster groups
    #cg.ax_col_dendrogram.clear()
    #for k, cat in enumerate(cluster_ids.categories):
    #    cg.ax_col_dendrogram.bar(0, 0, color=palettes[k], label=cat, linewidth=0)
    #cg.ax_col_dendrogram.legend(loc="center", ncol=15, fontsize=10)
    #cg.ax_col_dendrogram.grid(False)
    #cg.ax_col_dendrogram.set_xticks([])
    #cg.ax_col_dendrogram.set_yticks([])

    if not isinstance(data, anndata.AnnData):
        if cur_matkey != data.current_matrix():
            data.select_matrix(cur_matkey)

    return cg if not show else None


def dotplot(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    keys: Union[str, List[str]],
    by: str,
    reduce_function: Callable[[np.ndarray], float] = np.mean,
    fraction_min: float = 0,
    fraction_max: float = None,
    dot_min: int = 0,
    dot_max: int = 20,
    cmap: Union[str, List[str], Tuple[str]] = 'Reds',
    sort_function: Callable[[pd.DataFrame], List[str]] = None,
    grid: bool = True,
    show: bool = True,
    **kwds,
):
    """
    Generate a dot plot.

    Parameters
    ----------

    data: ``AnnData`` or ``UnimodalData`` or ``MultimodalData`` object
        Single cell expression data.
    keys: ``str`` or ``List[str]``
        Features to plot.
    by: ``str``
        Cell attribute to plot.
    reduce_function: ``Callable[[np.ndarray], float]``, optional, default: ``np.mean``
        Function to calculate statistic on expression data. Default is mean.
    fraction_min: ``float``, optional, default: ``0``.
        Minimum fraction of expressing cells to consider.
    fraction_max: ``float``, optional, default: ``None``.
        Maximum fraction of expressing cells to consider. If ``None``, use the maximum value from data.
    dot_min: ``int``, optional, default: ``0``.
        Minimum size in pixels for dots.
    dot_max: ``int``, optional, default: ``20``.
        Maximum size in pixels for dots.
    cmap: ``str`` or ``List[str]`` or ``Tuple[str]``, optional, default: ``Reds``
        Color map.
    sort_function: ``Callable[[pd.DataFrame], List[str]]``, optional, default: ``None``
        Function used for sorting labels. If ``None``, don't sort.
    grid: ``bool``, optional, default: ``True``
        If ``True``, plot grids.
    show: ``bool``, optional, default: ``True``
        Return a ``Figure`` object if ``False``; return ``None`` otherwise.
    **kwds:
        Are passed to ``matplotlib.pyplot.scatter``.

    Returns
    -------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``show == False``

    Examples
    --------
    >>> pg.dotplot(data, keys = ['CD14', 'TRAC', 'CD34'], by = 'louvain_labels')

    """
    sns.set(font_scale=0.7, style='whitegrid')

    if not is_list_like(keys):
        keys = [keys]

    keywords = dict(cmap=cmap)
    keywords.update(kwds)

    from scipy.sparse import issparse
    X = data[:, keys].X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(data=X, columns=keys)
    df[by] = data.obs[by].values

    def non_zero(g):
        return np.count_nonzero(g) / g.shape[0]

    summarized_df = df.groupby(by).aggregate([reduce_function, non_zero])
    if sort_function is not None:
        row_indices = sort_function(summarized_df)
        summarized_df = summarized_df.iloc[row_indices]
    else:
        summarized_df = summarized_df.loc[natsorted(summarized_df.index, reverse=True)]

    mean_columns = []
    frac_columns = []
    for j in range(len(summarized_df.columns)):
        if j % 2 == 0:
            mean_columns.append(summarized_df.columns[j])
        else:
            frac_columns.append(summarized_df.columns[j])

    # Features on rows, by on rows
    fraction_df = summarized_df[frac_columns]
    mean_df = summarized_df[mean_columns]

    y, x = np.indices(mean_df.shape)
    y = y.flatten()
    x = x.flatten()
    fraction = fraction_df.values.flatten()
    if fraction_max is None:
        fraction_max = fraction.max()
    pixels = _get_dot_size(fraction, fraction_min, fraction_max, dot_min, dot_max)
    summary_values = mean_df.values.flatten()
    xlabel = [keys[i] for i in range(len(keys))]
    ylabel = [str(summarized_df.index[i]) for i in range(len(summarized_df.index))]
    dotplot_df = pd.DataFrame(data=dict(x=x, y=y, value=summary_values, pixels=pixels, fraction=fraction,
            xlabel=np.array(xlabel)[x], ylabel=np.array(ylabel)[y]))

    xticks = keys
    yticks = summarized_df.index.map(str).values

    import matplotlib.gridspec as gridspec

    width = int(np.ceil(((dot_max + 1) + 4) * len(xticks) + dotplot_df['ylabel'].str.len().max()) + dot_max + 100)
    height = int(np.ceil(((dot_max + 1) + 4) * len(yticks) + dotplot_df['xlabel'].str.len().max()) + 50)
    fig = plt.figure(figsize=(1.1 * width / 100.0, height / 100.0), dpi=100)
    gs = gridspec.GridSpec(3, 11, figure = fig)

    # note we take the max label string length as an approximation of width of labels in pixels
    ax = fig.add_subplot(gs[:, :-1])

    sc = ax.scatter(x='x', y='y', c='value', s='pixels', data=dotplot_df, linewidth=0.5, edgecolors='black', **keywords)

    ax.spines["top"].set_color('black')
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.spines["right"].set_color('black')
    if not grid:
        ax.grid(False)

    ax.set_ylabel(str(by))
    ax.set_xlabel('')
    ax.set_xlim(-1, len(xticks))
    ax.set_ylim(-1, len(yticks))
    ax.set_xticks(range(len(xticks)))
    ax.set_xticklabels(xticks)
    ax.set_yticks(range(len(summarized_df.index)))
    ax.set_yticklabels(yticks)
    plt.xticks(rotation=90)

    cbar = plt.colorbar(sc)
    #cbar.set_label("Mean of\nexpressing cells")

    size_range = fraction_max - fraction_min
    if 0.3 < size_range <= 0.6:
        size_legend_step = 0.1
    elif size_range <= 0.3:
        size_legend_step = 0.05
    else:
        size_legend_step = 0.2

    size_ticks = np.arange(fraction_min if fraction_min > 0 or fraction_min > 0 else fraction_min + size_legend_step,
        fraction_max + size_legend_step, size_legend_step)

    ax2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0, -1])
    size_legend = fig.add_subplot(ax2[0])
    size_tick_pixels = _get_dot_size(size_ticks, fraction_min, fraction_max, dot_min, dot_max)

    size_tick_labels = ["{:.0%}".format(x) for x in size_ticks]
    size_legend.scatter(x=np.repeat(0, len(size_ticks)), y=np.arange(0, len(size_ticks)), s=size_tick_pixels, c='black', linewidth=0.5)
    size_legend.title.set_text("Fraction of\nexpressing cells")
    size_legend.set_xlim(-0.1, 0.1)
    size_legend.set_xticks([])

    ymin, ymax = size_legend.get_ylim()
    size_legend.set_ylim(ymin, ymax + 0.5)

    size_legend.set_yticks(np.arange(len(size_ticks)))
    size_legend.set_yticklabels(size_tick_labels)
    size_legend.tick_params(axis='y', labelleft=False, labelright=True)

    size_legend.spines["top"].set_visible(False)
    size_legend.spines["bottom"].set_visible(False)
    size_legend.spines["left"].set_visible(False)
    size_legend.spines["right"].set_visible(False)
    size_legend.grid(False)

    # Reset global settings.
    sns.set(font_scale=1.0, style='whitegrid')

    return fig if not show else None

def dendrogram(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    groupby: str,
    rep: str = 'pca',
    genes: Optional[List[str]] = None,
    correlation_method: str = 'pearson',
    n_clusters: Optional[int] = None,
    affinity: str = 'euclidean',
    linkage: str = 'complete',
    compute_full_tree: Union[str, bool] = 'auto',
    distance_threshold: Optional[float] = 0,
    figsize: Tuple[float, float] = (6, 6),
    orientation: str = 'top',
    color_threshold: Optional[float] = None,
    show: bool = True,
    **kwargs,
):
    """
    Generate a dendrogram on hierarchical clustering result.

    Parameters
    ----------

    data: ``MultimodalData``, ``UnimodalData``, or ``AnnData`` object
        Single cell expression data.
    groupby: ``str``
        Categorical cell attribute to plot, which must exist in ``data.obs``.
    genes: ``List[str]``, optional, default: ``None``
        List of genes to use. Gene names must exist in ``data.var``. If set, use the counts in ``data.X`` for plotting; if set as ``None``, use the embedding specified in ``rep`` for plotting.
    rep: ``str``, optional, default: ``pca``
        Cell embedding to use. It only works when ``genes``is ``None``, and its key ``"X_"+rep`` must exist in ``data.obsm``. By default, use PCA coordinates.
    correlation_method: ``str``, optional, default: ``pearson``
        Method of correlation between categories specified in ``data.obs``. Available options are: ``pearson``, ``kendall``, ``spearman``. See `pandas corr documentation`_ for details.
    n_clusters: ``int``, optional, default: ``None``
        The number of clusters to find, used by hierarchical clustering. It must be ``None`` if ``distance_threshold`` is not ``None``.
    affinity: ``str``, optional, default: ``correlation``
        Metric used to compute the linkage, used by hierarchical clustering. Valid values for metric are:
            - From scikit-learn: ``cityblock``, ``cosine``, ``euclidean``, ``l1``, ``l2``, ``manhattan``.
            - From scipy.spatial.distance: ``braycurtis``, ``canberra``, ``chebyshev``, ``correlation``, ``dice``, ``hamming``, ``jaccard``, ``kulsinski``, ``mahalanobis``, ``minkowski``, ``rogerstanimoto``, ``russellrao``, ``seuclidean``, ``sokalmichener``, ``sokalsneath``, ``sqeuclidean``, ``yule``.
        Default is the correlation distance. See `scikit-learn distance documentation <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html>`_ for details.
    linkage: ``str``, optional, default: ``complete``
        Which linkage criterion to use, used by hierarchical clustering. Below are available options:
            - ``ward`` minimizes the variance of the clusters being merged.
            - ``avarage`` uses the average of the distances of each observation of the two sets.
            - ``complete`` uses the maximum distances between all observations of the two sets. (Default)
            - ``single`` uses the minimum of the distances between all observations of the two sets.
        See `scikit-learn documentation`_ for details.
    compute_full_tree: ``str`` or ``bool``, optional, default: ``auto``
        Stop early the construction of the tree at ``n_clusters``, used by hierarchical clustering. It must be ``True`` if ``distance_threshold`` is not ``None``.
        By default, this option is ``auto``, which is ``True`` if and only if ``distance_threshold`` is not ``None``, or ``n_clusters`` is less than ``min(100, 0.02 * n_groups)``, where ``n_groups`` is the number of categories in ``data.obs[groupby]``.
    distance_threshold: ``float``, optional, default: ``0``
        The linkage distance threshold above which, clusters will not be merged. If not ``None``, ``n_clusters`` must be ``None`` and ``compute_full_tree`` must be ``True``.
    figsize: ``Tuple[float, float]``, optional, default: ``(6, 6)``
        The size (width, height) in inches of figure.
    orientation: ``str``, optional, default: ``top``
        The direction to plot the dendrogram. Available options are: ``top``, ``bottom``, ``left``, ``right``. See `scipy dendrogram documentation`_ for explanation.
    color_threshold: ``float``, optional, default: ``None``
        Threshold for coloring clusters. See `scipy dendrogram documentation`_ for explanation.
    show: ``bool``, optional, default: ``True``
        Return a ``Figure`` object if ``False``; return ``None`` otherwise.
    **kwargs:
        Are passed to ``scipy.cluster.hierarchy.dendrogram``.

    .. _scikit-learn documentation: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html
    .. _scipy dendrogram documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
    .. _pandas corr documentation: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.corr.html

    Returns
    -------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``show == False``

    Examples
    --------
    >>> pg.dendrogram(data, groupby='louvain_labels')
    >>> pg.dendrogram(data, groupby='louvain_labels', genes=data.var_names)
    """
    if genes is None:
        embed_df = pd.DataFrame(X_from_rep(data, rep))
        embed_df.set_index(data.obs[groupby], inplace=True)
    else:
        sub_data = data[:, genes]
        X = sub_data.X.toarray() if issparse(sub_data.X) else sub_data.X
        embed_df = pd.DataFrame(X)
        embed_df.set_index(data.obs[groupby], inplace=True)

    mean_df = embed_df.groupby(level=0).mean()

    from sklearn.cluster import AgglomerativeClustering
    from scipy.cluster.hierarchy import dendrogram

    corr_mat = mean_df.T.corr(method=correlation_method)

    clusterer = AgglomerativeClustering(
                    n_clusters=n_clusters,
                    affinity=affinity,
                    linkage=linkage,
                    compute_full_tree=compute_full_tree,
                    distance_threshold=distance_threshold
                )
    clusterer.fit(corr_mat)

    print(clusterer.distances_)

    counts = np.zeros(clusterer.children_.shape[0])
    n_samples = len(clusterer.labels_)
    for i, merge in enumerate(clusterer.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # Leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([clusterer.children_, clusterer.distances_, counts]).astype(float)

    fig, axis = plt.subplots(1, 1, figsize=figsize)
    dendrogram(linkage_matrix, labels=mean_df.index.categories, ax=axis, **kwargs)
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()

    return fig if not show else None
