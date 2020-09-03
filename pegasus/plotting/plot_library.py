import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.sparse import issparse
from pandas.api.types import is_numeric_dtype, is_categorical_dtype, is_list_like
from scipy.stats import zscore
from sklearn.metrics import adjusted_mutual_info_score
from natsort import natsorted


import anndata
from pegasusio import UnimodalData, MultimodalData

from typing import List, Tuple, Union, Optional, Callable

import logging
logger = logging.getLogger(__name__)

from pegasus.tools import X_from_rep
from .plot_utils import _transform_basis, _get_nrows_and_ncols, _get_marker_size, _get_dot_size, _get_subplot_layouts, _get_legend_ncol, _get_palette, _plot_cluster_labels_in_heatmap, RestrictionParser, DictWithDefault, _generate_categories, _plot_corners


def scatter(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attrs: Union[str, List[str]],
    basis: Optional[str] = "umap",
    matkey: Optional[str] = None,
    restrictions: Optional[Union[str, List[str]]] = None,
    show_background: Optional[bool] = False,
    alpha: Optional[Union[float, List[float]]] = 1.0,
    legend_loc: Optional[Union[str, List[str]]] = "right margin",
    legend_ncol: Optional[str] = None,
    palettes: Optional[Union[str, List[str]]] = None,
    cmaps: Optional[Union[str, List[str]]] = "YlOrRd",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    panel_size: Optional[Tuple[float, float]] = (4, 4),
    left: Optional[float] = 0.2,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.4,
    hspace: Optional[float] = 0.15,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
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
    restrictions: ``str`` or ``List[str]``, optional, default: None
        A list of restrictions to subset data for plotting. There are two types of restrictions: global restriction and attribute-specific restriction. Global restriction appiles to all attributes in ``attrs`` and takes the format of 'key:value,value...', or 'key:~value,value...'. This restriction selects cells with the ``data.obs[key]`` values belong to 'value,value...' (or not belong to if '~' shows). Attribute-specific restriction takes the format of 'attr:key:value,value...', or 'attr:key:~value,value...'. It only applies to one attribute 'attr'. If 'attr' and 'key' are the same, one can use '.' to replace 'key' (e.g. ``cluster_labels:.:value1,value2``).
    show_background: ``bool``, optional, default: False
        Only applicable if `restrictions` is set. By default, only data points selected are shown. If show_background is True, data points that are not selected will also be shown.
    alpha: ``float`` or ``List[float]``, optional, default: ``1.0``
        Alpha value for blending, from 0.0 (transparent) to 1.0 (opaque). If this is a list, the length must match attrs, which means we set a separate alpha value for each attribute.
    legend_loc: ``str`` or ``List[str]``, optional, default: ``right margin``
        Legend location. Can be either "right margin" or "on data". If a list is provided, set 'legend_loc' for each attribute in 'attrs' separately.
    legend_ncol: ``str``, optional, default: None
        Only applicable if legend_loc == "right margin". Set number of columns used to show legends.
    palettes: ``str`` or ``List[str]``, optional, default: None
        Used for setting colors for every categories in categorical attributes. Each string in ``palettes`` takes the format of 'attr:color1,color2,...,colorn'. 'attr' is the categorical attribute and 'color1' - 'colorn' are the colors for each category in 'attr' (e.g. 'cluster_labels:black,blue,red,...,yellow'). If there is only one categorical attribute in 'attrs', ``palletes`` can be set as a single string and the 'attr' keyword can be omitted (e.g. "blue,yellow,red").
    cmaps: ``str`` or ``List[str]``, optional, default: ``YlOrRd``
        Used for setting colormap for numeric attributes. Each string in ``cmaps`` takes the format of 'colormap' or 'attr:colormap'. 'colormap' sets the default colormap for all numeric attributes. 'attr:colormap' overwrites attribute 'attr's colormap as 'colormap'.
    vmin: ``float``, optional, default: None
        Minimum value to show on a numeric scatter plot (feature plot).
    vmax: ``float``, optional, default: None
        Maximum value to show on a numeric scatter plot (feature plot).
    nrows: ``int``, optional, default: None
        Number of rows in the figure. If not set, pegasus will figure it out automatically.
    ncols: ``int``, optional, default: None
        Number of columns in the figure. If not set, pegasus will figure it out automatically.
    panel_size: `tuple`, optional (default: `(6, 4)`)
        The panel size (width, height) in inches.
    left: `float`, optional (default: `0.2`)
        This parameter sets the figure's left margin as a fraction of panel's width (left * panel_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of panel's height (bottom * panel_size[1]).
    wspace: `float`, optional (default: `0.4`)
        This parameter sets the width between panels and also the figure's right margin as a fraction of panel's width (wspace * panel_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between panels and also the figure's top margin as a fraction of panel's height (hspace * panel_size[1]).
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: 300.0
        The resolution of the figure in dots-per-inch.

    Returns
    -------

    `Figure` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> pg.scatter(data, attrs=['louvain_labels', 'Channel'], basis='fitsne')
    >>> pg.scatter(data, attrs=['CD14', 'TRAC'], basis='umap')
    """
    if not is_list_like(attrs):
        attrs = [attrs]
    nattrs = len(attrs)

    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()

    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
        data.select_matrix(matkey)

    x = data.obsm[f"X_{basis}"][:, 0]
    y = data.obsm[f"X_{basis}"][:, 1]

    # four corners of the plot
    corners = np.array(np.meshgrid([x.min(), x.max()], [y.min(), y.max()])).T.reshape(-1, 2)


    basis = _transform_basis(basis)
    marker_size = _get_marker_size(x.size)
    nrows, ncols = _get_nrows_and_ncols(nattrs, nrows, ncols)
    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=ncols, panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=hspace, squeeze=False)


    if not is_list_like(alpha):
        alpha = [alpha] * nattrs

    if not is_list_like(legend_loc):
        legend_loc = [legend_loc] * nattrs
    legend_fontsize = [5 if x == 'on data' else 10 for x in legend_loc]

    palettes = DictWithDefault(palettes)
    cmaps = DictWithDefault(cmaps)
    restr_obj = RestrictionParser(restrictions)
    restr_obj.calc_default(data)

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])

            if i * ncols + j < nattrs:
                pos = i * ncols + j
                attr = attrs[pos]

                if attr in data.obs:
                    values = data.obs[attr].values
                else:
                    try:
                        loc = data.var_names.get_loc(attr)
                    except KeyError:
                        raise KeyError(f"{attr} is neither in data.obs nor data.var_names!")
                    values = data.X[:, loc].toarray().ravel() if issparse(data.X) else data.X[:, loc]

                selected = restr_obj.get_satisfied(data, attr)

                if is_numeric_dtype(values):
                    cmap = cmaps.get(attr, squeeze = True)
                    if cmap is None:
                        raise KeyError(f"Please set colormap for attribute {attr} or set a default colormap!")

                    _plot_corners(ax, corners, marker_size)

                    img = ax.scatter(
                        x[selected],
                        y[selected],
                        c=values[selected],
                        s=marker_size,
                        marker=".",
                        alpha=alpha[pos],
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
                    labels, with_background = _generate_categories(values, restr_obj.get_satisfied(data, attr))
                    label_size = labels.categories.size

                    palette = palettes.get(attr)
                    if palette is None:
                        palette = _get_palette(label_size, with_background=with_background, show_background=show_background)
                    elif with_background:
                        palette = ["gainsboro" if show_background else "white"] + list(palette)

                    text_list = []
                    for k, cat in enumerate(labels.categories):
                        idx = labels == cat
                        if idx.sum() > 0:
                            scatter_kwargs = {"marker": ".", "alpha": alpha[pos], "edgecolors": "none", "rasterized": True}

                            if cat != "":
                                if legend_loc[pos] != "on data":
                                    scatter_kwargs["label"] = cat
                                else:
                                    text_list.append((np.median(x[idx]), np.median(y[idx]), cat))

                            if cat != "" or (cat == "" and show_background):
                                ax.scatter(
                                    x[idx],
                                    y[idx],
                                    c=palette[k],
                                    s=marker_size,
                                    **scatter_kwargs,
                                )
                            else:
                                _plot_corners(ax, corners, marker_size)

                    if legend_loc[pos] == "right margin":
                        legend = ax.legend(
                            loc="center left",
                            bbox_to_anchor=(1, 0.5),
                            frameon=False,
                            fontsize=legend_fontsize[pos],
                            ncol=_get_legend_ncol(label_size, legend_ncol),
                        )
                        for handle in legend.legendHandles:
                            handle.set_sizes([300.0])
                    elif legend_loc[pos] == "on data":
                        texts = []
                        for px, py, txt in text_list:
                            texts.append(ax.text(px, py, txt, fontsize=legend_fontsize[pos], fontweight = "bold", ha = "center", va = "center"))
                        # from adjustText import adjust_text
                        # adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', lw=0.5))

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

    return fig if return_fig else None


def scatter_groups(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attr: str,
    groupby: str,
    basis: Optional[str] = "umap",
    matkey: Optional[str] = None,
    restrictions: Optional[Union[str, List[str]]] = None,
    show_full: Optional[bool] = True,
    categories: Optional[List[str]] = None,
    alpha: Optional[float] = 1.0,
    legend_loc: Optional[str] = "right margin",
    legend_ncol: Optional[str] = None,
    palette: Optional[str] = None,
    cmap: Optional[str] = "YlOrRd",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    nrows: Optional[int] = None,
    ncols: Optional[int] = None,
    panel_size: Optional[Tuple[float, float]] = (4, 4),
    left: Optional[float] = 0.2,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.4,
    hspace: Optional[float] = 0.15,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """ Generate scatter plots of attribute 'attr' for each category in attribute 'group'. Optionally show scatter plot containing data points from all categories in 'group'.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
       Use current selected modality in data.
    attr: ``str``
        Color scatter plots by attribute 'attr'. This attribute should be one key in data.obs or data.var_names (e.g. one gene). If it is categorical, a palette will be used to color each category separately. Otherwise, a color map will be used.
    groupby: ``str``
        Generate separate scatter plots of 'attr' for data points in each category in 'groupby', which should be a key in data.obs representing one categorical variable.
    basis: ``str``, optional, default: ``umap``
        Basis to be used to generate scatter plots. Can be either 'umap', 'tsne', 'fitsne', 'fle', 'net_tsne', 'net_fitsne', 'net_umap' or 'net_fle'.
    matkey: ``str``, optional, default: None
        If matkey is set, select matrix with matkey as keyword in the current modality. Only works for MultimodalData or UnimodalData objects.
    restrictions: ``str`` or ``List[str]``, optional, default: None
        A list of restrictions to subset data for plotting. Each restriction takes the format of 'key:value,value...', or 'key:~value,value...'. This restriction selects cells with the ``data.obs[key]`` values belong to 'value,value...' (or not belong to if '~' shows).
    show_full: ``bool``, optional, default: True
        Show the scatter plot with all categories in 'groupby' as the first plot.
    categories: ``List[str]``, optional, default: None
        Redefine group structure based on attribute 'groupby'. If 'categories' is not None, each string in the list takes the format of 'category_name:value,value', or 'category_name:~value,value...", where 'category_name' refers to new category name, 'value' refers to one of the category in 'groupby' and '~' refers to exclude values.
    alpha: ``float``, optional, default: ``1.0``
        Alpha value for blending, from 0.0 (transparent) to 1.0 (opaque).
    legend_loc: ``str``, optional, default: ``right margin``
        Legend location. Can be either "right margin" or "on data".
    legend_ncol: ``str``, optional, default: None
        Only applicable if legend_loc == "right margin". Set number of columns used to show legends.
    palette: ``str``, optional, default: None
        Used for setting colors for one categorical attribute (e.g. "black,blue,red,...,yellow").
    cmap: ``str``, optional, default: ``YlOrRd``
        Used for setting colormap for one numeric attribute.
    vmin: ``float``, optional, default: None
        Minimum value to show on a numeric scatter plot (feature plot).
    vmax: ``float``, optional, default: None
        Maximum value to show on a numeric scatter plot (feature plot).
    nrows: ``int``, optional, default: None
        Number of rows in the figure. If not set, pegasus will figure it out automatically.
    ncols: ``int``, optional, default: None
        Number of columns in the figure. If not set, pegasus will figure it out automatically.
    panel_size: `tuple`, optional (default: `(6, 4)`)
        The panel size (width, height) in inches.
    left: `float`, optional (default: `0.2`)
        This parameter sets the figure's left margin as a fraction of panel's width (left * panel_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of panel's height (bottom * panel_size[1]).
    wspace: `float`, optional (default: `0.4`)
        This parameter sets the width between panels and also the figure's right margin as a fraction of panel's width (wspace * panel_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between panels and also the figure's top margin as a fraction of panel's height (hspace * panel_size[1]).
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: 300.0
        The resolution of the figure in dots-per-inch.

    Returns
    -------

    `Figure` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> pg.scatter_groups(data, attr='louvain_labels', groupby='Individual', basis='tsne', nrows = 2, ncols = 4, alpha = 0.5)
    >>> pg.scatter_groups(data, attr='anno', groupby='Channel', basis='umap', categories=['new_cat1:channel1,channel2', 'new_cat2:channel3'])
    """
    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
        data.select_matrix(matkey)

    x = data.obsm[f"X_{basis}"][:, 0]
    y = data.obsm[f"X_{basis}"][:, 1]

    # four corners of the plot
    corners = np.array(np.meshgrid([x.min(), x.max()], [y.min(), y.max()])).T.reshape(-1, 2)

    basis = _transform_basis(basis)
    marker_size = _get_marker_size(x.size)

    if attr in data.obs:
        values = data.obs[attr].values
    else:
        try:
            loc = data.var_names.get_loc(attr)
        except KeyError:
            raise KeyError(f"{attr} is neither in data.obs nor data.var_names!")
        values = data.X[:, loc].toarray().ravel() if issparse(data.X) else data.X[:, loc]

    is_cat = is_categorical_dtype(values)
    if (not is_cat) and (not is_numeric_dtype(values)):
        values = pd.Categorical(values, categories=natsorted(np.unique(values)))
        is_cat = True

    assert groupby in data.obs
    groups = data.obs[groupby].values
    if not is_categorical_dtype(groups):
        groups = pd.Categorical(groups, categories=natsorted(np.unique(groups)))


    restr_obj = RestrictionParser(restrictions)
    restr_obj.calc_default(data)
    selected = restr_obj.get_satisfied(data)
    nsel = selected.sum()

    if nsel < data.shape[0]:
        x = x[selected]
        y = y[selected]
        values = values[selected]
        groups = groups[selected]

    df_g = pd.DataFrame()
    if show_full:
        df_g["All"] = np.ones(nsel, dtype=bool)
    if categories is None:
        for cat in groups.categories:
            df_g[cat] = groups == cat
    else:
        cat_obj = RestrictionParser(categories)
        for cat, idx in cat_obj.next_category(groups):
            df_g[cat] = idx

    nrows, ncols = _get_nrows_and_ncols(df_g.shape[1], nrows, ncols)
    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=ncols, panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=hspace, squeeze=False)

    legend_fontsize = 5 if legend_loc == 'on data' else 10

    if is_cat:
        labels = values
        label_size = labels.categories.size
        palette = _get_palette(label_size) if palette is None else np.array(palette.split(","))
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
                    text_list = []
                    for k, cat in enumerate(labels.categories):
                        idx = np.logical_and(df_g.iloc[:, gid].values, labels == cat)
                        _plot_corners(ax, corners, marker_size)
                        if idx.sum() > 0:
                            scatter_kwargs = {"marker": ".", "alpha": alpha, "edgecolors": "none", "rasterized": True}

                            if legend_loc != "on data":
                                scatter_kwargs["label"] = str(cat)
                            else:
                                text_list.append((np.median(x[idx]), np.median(y[idx]), str(cat)))

                            ax.scatter(
                                x[idx],
                                y[idx],
                                c=palette[k],
                                s=marker_size,
                                **scatter_kwargs,
                            )

                    if legend_loc == "right margin":
                        legend = ax.legend(
                            loc="center left",
                            bbox_to_anchor=(1, 0.5),
                            frameon=False,
                            fontsize=legend_fontsize,
                            ncol=legend_ncol,
                        )
                        for handle in legend.legendHandles:
                            handle.set_sizes([300.0])
                    elif legend_loc == "on data":
                        texts = []
                        for px, py, txt in text_list:
                            texts.append(ax.text(px, py, txt, fontsize=legend_fontsize, fontweight = "bold", ha = "center", va = "center"))
                else:
                    _plot_corners(ax, corners, marker_size)

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

    return fig if return_fig else None


def compo_plot(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    groupby: str,
    condition: str,
    style: Optional[str] = "frequency",
    restrictions: Optional[Union[str, List[str]]] = None,
    xlabel: Optional[str] = None,
    panel_size: Optional[Tuple[float, float]] = (6, 4),
    left: Optional[float] = 0.15,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.3,
    hspace: Optional[float] = 0.15,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """Generate a composition plot, which shows the percentage of cells from each condition for every cluster.

    This function is used to generate composition plots, which are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects.

    Parameters
    ----------

    data : `AnnData` or `UnimodalData` or `MultimodalData` object
        Single cell expression data.
    groupby : `str`
        A categorical variable in data.obs that is used to categorize the cells, e.g. cell type.
    condition: `str`
        A categorical variable in data.obs that is used to calculate frequency within each category defined by 'groupby', e.g. donor.
    style: `str`, optional (default: `frequency`)
        Composition plot style. Can be either `frequency`, or 'normalized'. Within each cluster, the `frequency` style show the percentage of cells from each 'condition' within each category in 'groupby' (stacked), the `normalized` style shows for each category in 'groupby' the percentage of cells that are also in each 'condition' over all cells that are in the same 'condition' (not stacked).
    restrictions: ``str`` or ``List[str]``, optional, default: None
        A list of restrictions to subset data for plotting. Each restriction takes the format of 'key:value,value...', or 'key:~value,value...'. This restriction selects cells with the ``data.obs[key]`` values belong to 'value,value...' (or not belong to if '~' shows).
    xlabel: `str`, optional (default None)
        Label for the horizontal axis. If None, use 'groupby'.
    panel_size: `tuple`, optional (default: `(6, 4)`)
        The plot size (width, height) in inches.
    left: `float`, optional (default: `0.15`)
        This parameter sets the figure's left margin as a fraction of panel's width (left * panel_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of panel's height (bottom * panel_size[1]).
    wspace: `float`, optional (default: `0.3`)
        This parameter sets the width between panels and also the figure's right margin as a fraction of panel's width (wspace * panel_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between panels and also the figure's top margin as a fraction of panel's height (hspace * panel_size[1]).
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    -------

    `Figure` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> fig = pg.compo_plot(data, 'louvain_labels', 'Donor', style = 'normalized')
    """
    if xlabel is None:
        xlabel = groupby

    fig, ax = _get_subplot_layouts(panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=hspace) # default nrows = 1 & ncols = 1

    restr_obj = RestrictionParser(restrictions)
    restr_obj.calc_default(data)
    selected = restr_obj.get_satisfied(data)

    df = pd.crosstab(data.obs.loc[selected, groupby], data.obs.loc[selected, condition])
    df = df.reindex(
        index=natsorted(df.index.values), columns=natsorted(df.columns.values)
    )

    if style == "frequency":
        df = df.div(df.sum(axis=1), axis=0) * 100.0
    else:
        assert style == "normalized"
        df = df.div(df.sum(axis=0), axis=1) * 100.0

    df.plot(
        kind = "bar",
        stacked = style == "frequency",
        legend = False,
        color = _get_palette(df.shape[1]),
        ax = ax,
    )

    ax.grid(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Percentage")
    ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5))

    if len(max(df.index.astype(str), key=len)) >= 5:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=-45, ha='left')

    return fig if return_fig else None


def violin(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    attrs: Union[str, List[str]],
    groupby: str,
    hue: Optional[str] = None,
    matkey: Optional[str] = None,
    stripplot: bool = False,
    scale: str = 'width',
    panel_size: Optional[Tuple[float, float]] = (8, 0.5),
    left: Optional[float] = 0.15,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.1,
    ylabel: Optional[str] = None,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """
    Generate a stacked violin plot.

    Parameters
    ----------
    data: ``AnnData`` or ``MultimodalData`` or ``UnimodalData`` object
        Single-cell expression data.
    attrs: ``str`` or ``List[str]``
        Cell attributes or features to plot.
        Cell attributes must exist in ``data.obs`` and must be numeric.
        Features must exist in ``data.var``.
    groupby: ``str``
        A categorical variable in data.obs that is used to categorize the cells, e.g. Clusters.
    hue: ``str``, optional, default: None
        'hue' should be a categorical variable in data.obs that has only two levels. Set 'hue' will show us split violin plots.
    matkey: ``str``, optional, default: ``None``
        If matkey is set, select matrix with matkey as keyword in the current modality. Only works for MultimodalData or UnimodalData objects.
    stripplot: ``bool``, optional, default: ``False``
        Attach a stripplot to the violinplot or not. This option will be automatically turn off if 'hue' is set.
    scale: ``str``, optional, default: ``width``
        The method used to scale the width of each violin:
            - If ``width``, each violin will have the same width.
            - If ``area``, each violin will have the same area.
            - If ``count``, the width of the violins will be scaled by the number of observations in that bin.
    jitter: ``float`` or ``bool``, optional, default: ``False``
        Amount of jitter (only along the categorical axis) to apply to stripplot. This is used only when ``stripplot`` is set to ``True``.
        This can be useful when you have many points and they overlap, so that it is easier to see the distribution. You can specify the amount of jitter (half the width of the uniform random variable support), or just use ``True`` for a good default.
    panel_size: ``Tuple[float, float]``, optional, default: ``(8, 0.5)``
        The size (width, height) in inches of each violin panel.
    left: ``float``, optional, default: ``0.15``
        This parameter sets the figure's left margin as a fraction of panel's width (left * panel_size[0]).
    bottom: ``float``, optional, default: ``0.15``
        This parameter sets the figure's bottom margin as a fraction of panel's height (bottom * panel_size[1]).
    wspace: ``float``, optional, default: ``0.1``
        This parameter sets the width between panels and also the figure's right margin as a fraction of panel's width (wspace * panel_size[0]).
    ylabel: ``str``, optional, default: ``None``
        Y-axis label. No label to show if ``None``.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.
    kwargs
        Are passed to ``seaborn.violinplot``.

    Returns
    -------
    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``show == False``

    Examples
    --------
    >>> pg.violin(data, attrs=['CD14', 'TRAC', 'CD34'], groupby='louvain_labels')
    """
    if not is_list_like(attrs):
        attrs = [attrs]

    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
        data.select_matrix(matkey)

    nrows = len(attrs)
    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=1, panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=0, squeeze=False, sharey=False)

    obs_keys = []
    genes = []

    for key in attrs:
        if key in data.obs:
            assert is_numeric_dtype(data.obs[key])
            obs_keys.append(key)
        else:
            if key not in data.var_names:
                logger.warning("Cannot find gene {key}. Please make sure all genes are included in data.var_names before running this function!")
                return None
            genes.append(key)

    df_list = [pd.DataFrame({"label": data.obs[groupby].values})]
    if hue is not None:
        df_list.append(pd.DataFrame({hue: data.obs[hue].values}))
        stripplot = False
    if len(obs_keys) > 0:
        df_list.append(data.obs[obs_keys].reset_index(drop=True))
    if len(genes) > 0:
        expr_mat = data[:, genes].X.toarray()
        df_list.append(pd.DataFrame(data=expr_mat, columns=genes))
    df = pd.concat(df_list, axis = 1)

    for i in range(nrows):
        ax = axes[i, 0]
        if stripplot:
            sns.stripplot(x="label", y=attrs[i], hue = hue, data=df, ax=ax, size=1, color="k", jitter=True)
        sns.violinplot(x="label", y=attrs[i], hue = hue, data=df, inner=None, linewidth=1, ax=ax, cut=0, scale=scale, split=True, **kwargs)
        ax.grid(False)

        if hue is not None:
            if i == 0:
                ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
            else:
                ax.get_legend().set_visible(False)

        if i < nrows - 1:
            ax.set_xlabel("")
        else:
            ax.set_xlabel(groupby)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_ylabel(attrs[i], labelpad=8, rotation=0, horizontalalignment='right', fontsize='medium')
        ax.tick_params(axis='y', right=True, left=False, labelright=True, labelleft=False, labelsize='small')

    if ylabel is not None:
        fig.text(0.02, 0.5, ylabel, rotation="vertical", fontsize="xx-large")

    # Reset current matrix if needed.
    if not isinstance(data, anndata.AnnData):
        if data.current_matrix() != cur_matkey:
            data.select_matrix(cur_matkey)

    return fig if return_fig else None


def heatmap(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    genes: Union[str, List[str]],
    groupby: str,
    matkey: Optional[str] = None,
    on_average: bool = True,
    switch_axes: bool = False,
    row_cluster: Optional[bool] = None,
    col_cluster: Optional[bool] = None,
    panel_size: Tuple[float, float] = (10, 10),
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """
    Generate a heatmap.

    Parameters
    -----------

    data: ``AnnData`` or ``MultimodalData`` or ``UnimodalData`` object
        Single-cell expression data.
    genes: ``str`` or ``List[str]``
        Features to plot.
    groupby: ``str``
        A categorical variable in data.obs that is used to categorize the cells, e.g. Clusters.
    matkey: ``str``, optional, default: ``None``
        If matkey is set, select matrix with matkey as keyword in the current modality. Only works for MultimodalData or UnimodalData objects.
    on_average: ``bool``, optional, default: ``True``
        If ``True``, plot cluster average gene expression (i.e. show a Matrixplot); otherwise, plot a general heatmap.
    switch_axes: ``bool``, optional, default: ``False``
        By default, X axis is for genes, and Y axis for clusters. If this parameter is ``True``, switch the axes.
        Moreover, with ``on_average`` being ``False``, if ``switch_axes`` is ``False``, ``row_cluster`` is enforced to be ``False``; if ``switch_axes`` is ``True``, ``col_cluster`` is enforced to be ``False``.
    row_cluster: ``bool``, optional, default: ``False``
        Cluster rows and generate a row-wise dendrogram.
    col_cluster: ``bool``, optional, default: ``True``
        Cluster columns and generate a column-wise dendrogram.
    panel_size: ``Tuple[float, float]``, optional, default: ``(10, 10)``
        Overall size of the heatmap in ``(width, height)`` form.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.
    kwargs
        Are passed to ``seaborn.heatmap``.

    .. _colormap documentation: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    Returns
    -------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> pg.heatmap(data, genes=['CD14', 'TRAC', 'CD34'], groupby='louvain_labels')

    """
    if isinstance(data, MultimodalData) or isinstance(data, UnimodalData):
        cur_matkey = data.current_matrix()
    if matkey is not None:
        assert isinstance(data, MultimodalData) or isinstance(data, UnimodalData)
        data.select_matrix(matkey)

    if row_cluster is None:
        row_cluster = True if switch_axes else False

    if col_cluster is None:
        col_cluster = True if not switch_axes else False

    df = pd.DataFrame(data[:, genes].X.toarray(), index=data.obs.index, columns=genes)
    df['cluster_name'] = data.obs[groupby]

    if on_average:
        if not 'cmap' in kwargs.keys():
            kwargs['cmap'] = 'Reds'
        df = df.groupby('cluster_name').mean()
        cluster_ids = df.index
    else:
        row_cluster = False if not switch_axes else row_cluster
        col_cluster = False if switch_axes else col_cluster

        cluster_ids = pd.Categorical(data.obs[groupby])
        idx = cluster_ids.argsort()
        df = df.iloc[idx, :]  # organize df by category order
        df.drop(columns=['cluster_name'], inplace=True)

        cell_colors = np.zeros(df.shape[0], dtype=object)
        palette = _get_palette(cluster_ids.categories.size)
        cluster_ids = cluster_ids[idx]
        for k, cat in enumerate(cluster_ids.categories):
            cell_colors[np.isin(cluster_ids, cat)] = palette[k]

    if not switch_axes:
        cg = sns.clustermap(
            data=df,
            row_colors=cell_colors if not on_average else None,
            col_colors=None,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            linewidths=0,
            yticklabels=cluster_ids if on_average else [],
            xticklabels=genes,
            figsize=panel_size,
            **kwargs,
        )
        cg.ax_heatmap.set_ylabel("")
    else:
        cg = sns.clustermap(
            data=df.T,
            row_colors=None,
            col_colors=cell_colors if not on_average else None,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            linewidths=0,
            yticklabels=genes,
            xticklabels=cluster_ids if on_average else [],
            figsize=panel_size,
            **kwargs,
        )
        cg.ax_heatmap.set_xlabel("")

    if row_cluster:
        cg.ax_heatmap.yaxis.tick_right()
    else:
        cg.ax_heatmap.yaxis.tick_left()

    cg.ax_row_dendrogram.set_visible(row_cluster)
    cg.cax.tick_params(labelsize=10)
    cg.fig.dpi = dpi

    if not row_cluster:
        # Move the colorbar to the right-side.
        color_box = cg.ax_heatmap.get_position()
        color_box.x0 = color_box.x1 + 0.04
        color_box.x1 = color_box.x0 + 0.02
        cg.cax.set_position(color_box)
        cg.cax.yaxis.set_ticks_position("right")
    else:
        # Avoid overlap of colorbar and row dendrogram.
        color_box = cg.cax.get_position()
        square_plot = cg.ax_heatmap.get_position()
        if square_plot.y1 > color_box.y0:
            y_diff = square_plot.y1 - color_box.y0
            color_box.y0 = square_plot.y1
            color_box.y1 += y_diff
            cg.cax.set_position(color_box)

    if not on_average:
        orientation = 'left' if not switch_axes else 'top'
        _plot_cluster_labels_in_heatmap(cg.ax_heatmap, cluster_ids, orientation)

    if not isinstance(data, anndata.AnnData):
        if cur_matkey != data.current_matrix():
            data.select_matrix(cur_matkey)

    return cg.fig if return_fig else None


def dotplot(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    genes: Union[str, List[str]],
    groupby: str,
    reduce_function: Callable[[np.ndarray], float] = np.mean,
    fraction_min: float = 0,
    fraction_max: float = None,
    dot_min: int = 0,
    dot_max: int = 20,
    switch_axes: bool = False,
    cmap: Union[str, List[str], Tuple[str]] = 'Reds',
    sort_function: Callable[[pd.DataFrame], List[str]] = None,
    grid: bool = True,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwds,
) -> Union[plt.Figure, None]:
    """
    Generate a dot plot.

    Parameters
    ----------

    data: ``AnnData`` or ``UnimodalData`` or ``MultimodalData`` object
        Single cell expression data.
    genes: ``str`` or ``List[str]``
        Features to plot.
    groupby: ``str``
        A categorical variable in data.obs that is used to categorize the cells, e.g. Clusters.
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
    switch_axes: ``bool``, optional, default: ``False``.
        If ``True``, switch X and Y axes.
    cmap: ``str`` or ``List[str]`` or ``Tuple[str]``, optional, default: ``Reds``
        Color map.
    sort_function: ``Callable[[pd.DataFrame], List[str]]``, optional, default: ``None``
        Function used for sorting labels. If ``None``, don't sort.
    grid: ``bool``, optional, default: ``True``
        If ``True``, plot grids.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.
    **kwds:
        Are passed to ``matplotlib.pyplot.scatter``.

    Returns
    -------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> pg.dotplot(data, genes = ['CD14', 'TRAC', 'CD34'], groupby = 'louvain_labels')

    """
    sns.set(font_scale=0.7, style='whitegrid')

    if not is_list_like(genes):
        geness = [genes]

    keywords = dict(cmap=cmap)
    keywords.update(kwds)

    from scipy.sparse import issparse
    X = data[:, genes].X
    if issparse(X):
        X = X.toarray()
    df = pd.DataFrame(data=X, columns=genes)
    df[groupby] = data.obs[groupby].values

    def non_zero(g):
        return np.count_nonzero(g) / g.shape[0]

    summarized_df = df.groupby(groupby).aggregate([reduce_function, non_zero])
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

    # Genes on columns, groupby on rows
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

    xlabel = [genes[i] for i in range(len(genes))]
    ylabel = [str(summarized_df.index[i]) for i in range(len(summarized_df.index))]

    xticks = genes
    yticks = summarized_df.index.map(str).values

    if switch_axes:
        x, y = y, x
        xlabel, ylabel = ylabel, xlabel
        xticks, yticks = yticks, xticks

    dotplot_df = pd.DataFrame(data=dict(x=x, y=y, value=summary_values, pixels=pixels, fraction=fraction,
                    xlabel=np.array(xlabel)[x], ylabel=np.array(ylabel)[y]))

    import matplotlib.gridspec as gridspec

    width = int(np.ceil(((dot_max + 1) + 4) * len(xticks) + dotplot_df['ylabel'].str.len().max()) + dot_max + 100)
    height = int(np.ceil(((dot_max + 1) + 4) * len(yticks) + dotplot_df['xlabel'].str.len().max()) + 50)
    fig = plt.figure(figsize=(1.1 * width / 100.0, height / 100.0), dpi=dpi)
    gs = gridspec.GridSpec(3, 11, figure = fig)

    # Main plot
    mainplot_col_grid = -2 if len(xlabel) < 10 else -1
    ax = fig.add_subplot(gs[:, :mainplot_col_grid])

    sc = ax.scatter(x='x', y='y', c='value', s='pixels', data=dotplot_df, linewidth=0.5, edgecolors='black', **keywords)

    ax.spines["top"].set_color('black')
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.spines["right"].set_color('black')
    if not grid:
        ax.grid(False)

    if not switch_axes:
        ax.set_ylabel(str(groupby))
        ax.set_xlabel('')
    else:
        ax.set_ylabel('')
        ax.set_xlabel(str(groupby))

    ax.set_xlim(-1, len(xticks))
    ax.set_ylim(-1, len(yticks))
    ax.set_xticks(range(len(xticks)))
    ax.set_xticklabels(xticks)
    ax.set_yticks(range(len(yticks)))
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

    legend_row_grid = 1 if height / 3 > 100 else 3
    ax2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0:legend_row_grid, -1])
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
    sns.reset_orig()

    return fig if return_fig else None


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
    panel_size: Tuple[float, float] = (6, 6),
    orientation: str = 'top',
    color_threshold: Optional[float] = None,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """
    Generate a dendrogram on hierarchical clustering result.

    The metrics used here are consistent with SCANPY's dendrogram_ implementation.

    *scikit-learn* `Agglomerative Clustering`_ implementation is used for hierarchical clustering.

    .. _dendrogram: https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dendrogram.html
    .. _Agglomerative Clustering: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html

    Parameters
    ----------

    data: ``MultimodalData``, ``UnimodalData``, or ``AnnData`` object
        Single cell expression data.
    genes: ``List[str]``, optional, default: ``None``
        List of genes to use. Gene names must exist in ``data.var``. If set, use the counts in ``data.X`` for plotting; if set as ``None``, use the embedding specified in ``rep`` for plotting.
    rep: ``str``, optional, default: ``pca``
        Cell embedding to use. It only works when ``genes``is ``None``, and its key ``"X_"+rep`` must exist in ``data.obsm``. By default, use PCA coordinates.
    groupby: ``str``
        Categorical cell attribute to plot, which must exist in ``data.obs``.
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
    panel_size: ``Tuple[float, float]``, optional, default: ``(6, 6)``
        The size (width, height) in inches of figure.
    orientation: ``str``, optional, default: ``top``
        The direction to plot the dendrogram. Available options are: ``top``, ``bottom``, ``left``, ``right``. See `scipy dendrogram documentation`_ for explanation.
    color_threshold: ``float``, optional, default: ``None``
        Threshold for coloring clusters. See `scipy dendrogram documentation`_ for explanation.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.
    **kwargs:
        Are passed to ``scipy.cluster.hierarchy.dendrogram``.

    .. _scikit-learn documentation: https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html
    .. _scipy dendrogram documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
    .. _pandas corr documentation: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.corr.html

    Returns
    -------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> pg.dendrogram(data, genes=data.var_names, groupby='louvain_labels')
    >>> pg.dendrogram(data, rep='pca', groupby='louvain_labels')
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
    mean_df.index = mean_df.index.astype('category')

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

    fig, ax = _get_subplot_layouts(panel_size=panel_size, dpi=dpi)
    dendrogram(linkage_matrix, labels=mean_df.index.categories, ax=ax, **kwargs)
    plt.xticks(rotation=90, fontsize=10)
    plt.tight_layout()

    return fig if return_fig else None


def hvfplot(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    top_n: int = 20,
    panel_size: Optional[Tuple[float, float]] = (6, 4),
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
) -> Union[plt.Figure, None]:
    """
    Generate highly variable feature plot.
    Only works for HVGs returned by ``highly_variable_features`` method with ``flavor=='pegasus'``.

    Parameters
    -----------

    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData`` object.
        Single cell expression data.
    top_n: ``int``, optional, default: ``20``
        Number of top highly variable features to show names.
    panel_size: ``Tuple[float, float]``, optional, default: ``(6, 4)``
        The size (width, height) in inches of figure.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    --------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    ---------
    >>> pg.hvfplot(data)
    >>> pg.hvfplot(data, top_n=10, dpi=150)
    """
    robust_idx = data.var["robust"].values
    x = data.var.loc[robust_idx, "mean"]
    y = data.var.loc[robust_idx, "var"]
    fitted = data.var.loc[robust_idx, "hvf_loess"]
    hvg_index = data.var.loc[robust_idx, "highly_variable_features"]
    hvg_rank = data.var.loc[robust_idx, "hvf_rank"]
    gene_symbols = data.var_names[robust_idx]

    fig, ax = _get_subplot_layouts(panel_size=panel_size, dpi=dpi)

    ax.scatter(x[hvg_index], y[hvg_index], s=5, c='b', marker='o', linewidth=0.5, alpha=0.5, label='highly variable features')
    ax.scatter(x[~hvg_index], y[~hvg_index], s=5, c='k', marker='o', linewidth=0.5, alpha=0.5, label = 'other features')
    ax.legend(loc = 'upper right', fontsize = 5)
    ax.set_xlabel("Mean log expression")
    ax.set_ylabel("Variance of log expression")

    order = x.argsort().values
    ax.plot(x[order], fitted[order], "r-", linewidth=1)

    ord_rank = hvg_rank.argsort().values

    texts = []
    for i in range(top_n):
        pos = ord_rank[i]
        texts.append(ax.text(x[pos], y[pos], gene_symbols[pos], fontsize=5))

    from adjustText import adjust_text
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', lw=0.5))

    return fig if return_fig else None


def qcviolin(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    plot_type: str,
    min_genes_before_filt: Optional[int] = 100,
    n_violin_per_panel: Optional[int] = 8,
    panel_size: Optional[Tuple[float, float]] = (6, 4),
    left: Optional[float] = 0.2,
    bottom: Optional[float] = 0.15,
    wspace: Optional[float] = 0.3,
    hspace: Optional[float] = 0.35,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
) -> Union[plt.Figure, None]:
    """
    Plot quality control statistics (before filtration vs. after filtration) as violin plots. Require statistics such as "n_genes", "n_counts" and "percent_mito" precomputed.

    Parameters
    -----------

    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData`` object.
        Single cell expression data.
    plot_type: ``str``
        Choose from ``gene``, ``count`` and ``mito``, which shows number of expressed genes, number of UMIs and percentage of mitochondrial rate.
    min_genes_before_filt: ``int``, optional, default: 100
        If data loaded are raw data (i.e. min(n_genes) == 0), filter out cell barcodes with less than ``min_genes_before_filt`` for better visual effects.
    n_violin_per_panel: ``int``, optional, default: 8
        Number of violin plots (samples) shown in one panel.
    panel_size: `tuple`, optional (default: `(6, 4)`)
        The panel size (width, height) in inches.
    left: `float`, optional (default: `0.2`)
        This parameter sets the figure's left margin as a fraction of panel's width (left * panel_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of panel's height (bottom * panel_size[1]).
    wspace: `float`, optional (default: `0.4`)
        This parameter sets the width between panels and also the figure's right margin as a fraction of panel's width (wspace * panel_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between panels and also the figure's top margin as a fraction of panel's height (hspace * panel_size[1]).
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    --------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    ---------
    >>> pg.qcviolin(data, "mito", dpi = 500)
    """
    pt2attr = {"gene": "n_genes", "count": "n_counts", "mito": "percent_mito"}
    pt2ylab = {
        "gene": "Number of expressed genes",
        "count": "Number of UMIs",
        "mito": "Percentage of mitochondrial UMIs",
    }

    if "df_qcplot" not in data.uns:
        if "Channel" not in data.obs:
            data.obs["Channel"] = pd.Categorical([""] * data.shape[0])

        target_cols = np.array(["Channel", "n_genes", "n_counts", "percent_mito"])
        target_cols = target_cols[np.isin(target_cols, data.obs.columns)]

        df = data.obs[data.obs["n_genes"] >= min_genes_before_filt] if data.obs["n_genes"].min() == 0 else data.obs
        df_plot_before = df[target_cols].copy()
        df_plot_before.reset_index(drop=True, inplace=True)
        df_plot_before["status"] = "original"

        df_plot_after = data.obs.loc[data.obs["passed_qc"], target_cols].copy()
        df_plot_after.reset_index(drop=True, inplace=True)
        df_plot_after["status"] = "filtered"

        df_qcplot = pd.concat((df_plot_before, df_plot_after), axis=0)

        df_qcplot["status"] = pd.Categorical(df_qcplot["status"].values, categories = ["original", "filtered"])
        df_qcplot["Channel"] = pd.Categorical(df_qcplot["Channel"].values, categories = natsorted(df_qcplot["Channel"].astype(str).unique()))

        data.uns["df_qcplot"] = df_qcplot


    df_qcplot = data.uns["df_qcplot"]

    if pt2attr[plot_type] not in df_qcplot:
        logger.warning(f"Cannot find qc metric {pt2attr[plot_type]}!")
        return None

    channels = df_qcplot["Channel"].cat.categories
    n_channels = channels.size
    n_pannels = (n_channels - 1) // n_violin_per_panel + 1

    nrows = ncols = None
    nrows, ncols = _get_nrows_and_ncols(n_pannels, nrows, ncols)
    fig, axes = _get_subplot_layouts(nrows=nrows, ncols=ncols, panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=hspace, sharex = False, sharey = False, squeeze=False)

    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            ax.grid(False)
            panel_no = i * ncols + j
            if panel_no < n_pannels:
                start = panel_no * n_violin_per_panel
                end = min(start + n_violin_per_panel, n_channels)
                idx = np.isin(df_qcplot["Channel"], channels[start:end])

                if start == 0 and end == n_channels:
                    df_plot = df_qcplot
                else:
                    df_plot = df_qcplot[idx].copy()
                    df_plot["Channel"] = pd.Categorical(df_plot["Channel"].values, categories = natsorted(channels[start:end]))

                sns.violinplot(
                    x="Channel",
                    y=pt2attr[plot_type],
                    hue="status",
                    data=df_plot,
                    split=True,
                    linewidth=0.5,
                    cut=0,
                    inner=None,
                    ax = ax,
                )

                ax.set_xlabel("Channel")
                ax.set_ylabel(pt2ylab[plot_type])
                ax.legend(loc="upper right", fontsize=8)
                if max([len(x) for x in channels[start:end]]) >= 5:
                    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=-45)
            else:
                ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])

    return fig if return_fig else None


def volcano(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    cluster_id: str,
    de_key: str = "de_res",
    de_test: str = 'mwu',
    qval_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    top_n: int = 20,
    panel_size: Optional[Tuple[float, float]] = (6, 4),
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
) -> Union[plt.Figure, None]:
    """
    Generate Volcano plots (-log10 p value vs. log2 fold change) for visualizing DE results.

    Parameters
    -----------

    data: ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData`` object.
        Single cell expression data.
    cluster_id: ``str``
        Cluster ID for the cluster we want to show DE results. There are two cases:
            * If ``condition`` is ``None`` in ``pg.de_analysis``: Just specify one cluster
              label in the cluster attribute used in ``pg.de_analysis``.
            * If ``condition`` is not ``None`` in ``pg.de_analysis``: Specify cluster ID in
              this format: **"cluster_label:cond_level"**, where **cluster_label** is the
              cluster label, and **cond_level** is the condition ID. And this shows result of
              cells within the cluster under the specific condition.
    de_key: ``str``, optional, default: ``de_res``
        The varm keyword for DE results. data.varm[de_key] should store the full DE result table.
    de_test: ``str``, optional, default: ``mwu``
        Which DE test results to show. Use MWU test result by default.
    qval_threshold: ``float``, optional, default: 0.05.
        Selected FDR rate. A horizontal line indicating this rate will be shown in the figure.
    log2fc_threshold: ``float``, optional, default: 1.0
        Log2 fold change threshold to highlight biologically interesting genes. Two vertical lines representing negative and positive log2 fold change will be shown.
    top_n: ``int``, optional, default: ``20``
        Number of top DE genes to show names. Genes are ranked by Log2 fold change.
    panel_size: ``Tuple[float, float]``, optional, default: ``(6, 4)``
        The size (width, height) in inches of figure.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    --------

    ``Figure`` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    ---------
    >>> pg.volcano(data, cluster_id = '1', dpi=200)
    """
    if de_key not in data.varm:
        logger.warning(f"Cannot find DE results '{de_key}'. Please conduct DE analysis first!")
        return None

    de_res = data.varm[de_key]

    fcstr = f"{cluster_id}:log2FC"
    pstr = f"{cluster_id}:{de_test}_pval"
    qstr = f"{cluster_id}:{de_test}_qval"

    columns = de_res.dtype.names
    if (fcstr not in columns) or (pstr not in columns) or (qstr not in columns):
        logger.warning(f"Please conduct DE test {de_test} first!")
        return None

    log2fc = de_res[fcstr]
    pvals = de_res[pstr]
    pvals[pvals == 0.0] = 1e-45 # very small pvalue to avoid log10 0
    neglog10p = -np.log10(pvals)
    yconst = min(neglog10p[de_res[qstr] <= qval_threshold])

    fig, ax = _get_subplot_layouts(panel_size=panel_size, dpi=dpi)

    idxsig = neglog10p >= yconst
    idxnsig = neglog10p < yconst
    idxfc = (log2fc <= -log2fc_threshold) | (log2fc >= log2fc_threshold)
    idxnfc = ~idxfc

    idx = idxnsig & idxnfc
    ax.scatter(log2fc[idx], neglog10p[idx], s=5, c='k', marker='o', linewidths=0.5, alpha=0.5, label="NS")
    idx = idxnsig & idxfc
    ax.scatter(log2fc[idx], neglog10p[idx], s=5, c='g', marker='o', linewidths=0.5, alpha=0.5, label=r"Log$_2$ FC")
    idx = idxsig & idxnfc
    ax.scatter(log2fc[idx], neglog10p[idx], s=5, c='b', marker='o', linewidths=0.5, alpha=0.5, label=r"q-value")
    idx = idxsig & idxfc
    ax.scatter(log2fc[idx], neglog10p[idx], s=5, c='r', marker='o', linewidths=0.5, alpha=0.5, label=r"q-value and log$_2$ FC")

    ax.set_xlabel(r"Log$_2$ fold change")
    ax.set_ylabel(r"$-$Log$_{10}$ $P$")

    legend = ax.legend(
        loc="center",
        bbox_to_anchor=(0.5, 1.1),
        frameon=False,
        fontsize=8,
        ncol=4,
    )
    for handle in legend.legendHandles: # adjust legend size
        handle.set_sizes([50.0])

    ax.axhline(y = yconst, c = 'k', lw = 0.5, ls = '--')
    ax.axvline(x = -log2fc_threshold, c = 'k', lw = 0.5, ls = '--')
    ax.axvline(x = log2fc_threshold, c = 'k', lw = 0.5, ls = '--')

    texts = []

    idx = np.where(idxsig & (log2fc >= log2fc_threshold))[0]
    posvec = np.argsort(log2fc[idx])[::-1][0:top_n]
    for pos in posvec:
        gid = idx[pos]
        texts.append(ax.text(log2fc[gid], neglog10p[gid], data.var_names[gid], fontsize=5))

    idx = np.where(idxsig & (log2fc <= -log2fc_threshold))[0]
    posvec = np.argsort(log2fc[idx])[0:top_n]
    for pos in posvec:
        gid = idx[pos]
        texts.append(ax.text(log2fc[gid], neglog10p[gid], data.var_names[gid], fontsize=5))

    from adjustText import adjust_text
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', lw=0.5))

    return fig if return_fig else None


def rank_plot(
    data: Union[MultimodalData, UnimodalData, anndata.AnnData],
    panel_size: Optional[Tuple[float, float]] = (6, 4),
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """Generate a barcode rank plot, which shows the total UMIs against barcode rank (in descending order with respect to total UMIs)

    Parameters
    ----------

    data : `AnnData` or `UnimodalData` or `MultimodalData` object
        The main data object.
    panel_size: `tuple`, optional (default: `(6, 4)`)
        The plot size (width, height) in inches.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    -------

    `Figure` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> fig = pg.rank_plot(data, show = False, dpi = 500)
    """
    fig, ax = _get_subplot_layouts(panel_size=panel_size, dpi=dpi) # default nrows = 1 & ncols = 1

    numis = data.X.sum(axis = 1).A1
    ords = np.argsort(numis)[::-1]
    ranks = np.array(range(1, numis.size + 1))
    ax.scatter(ranks, numis[ords], c = 'lightgrey', s = 5)
    ax.set_xscale("log", basex = 10)
    ax.set_yscale("log", basey = 10)
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("Total UMIs")

    def _gen_ticklabels(ticks, max_value):
        label_arr = ['1', '10', '100', '1000', '10K', '100K', '1M', '10M', '100M']
        ticklabels = [''] * ticks.size
        for i in range(ticks.size):
            exponent = int(round(np.log10(ticks[i])))
            if exponent >= 0 and ticks[i] <= max_value:
                ticklabels[i] = label_arr[exponent]
        return ticklabels

    ax.set_xticklabels(_gen_ticklabels(ax.get_xticks(), numis.size))
    ax.set_yticklabels(_gen_ticklabels(ax.get_yticks(), numis.max()))

    return fig if return_fig else None


def ridgeplot(
    data: Union[MultimodalData, UnimodalData],
    features: Union[str, List[str]],
    qc_attr: Optional[str] = None,
    overlap: Optional[float] = 0.5,
    panel_size: Optional[Tuple[float, float]] = (6, 4),
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """Generate ridge plots, up to 8 features can be shown in one figure.

    Parameters
    ----------

    data : `UnimodalData` or `MultimodalData` object
        CITE-Seq or Cyto data.
    features : `str` or `List[str]`
        One or more features to display.
    qc_attr: `str`, optional, default None
        If not None, only data.obs[qc_attr] == True are used.
    overlap: `float`, default 0.5
        Overlap between adjacent ridge plots (top and bottom).
    panel_size: `tuple`, optional (default: `(6, 4)`)
        The plot size (width, height) in inches.
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    -------

    `Figure` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    >>> fig = pg.ridgeplot(data, features = ['CD8', 'CD4', 'CD3'], show = False, dpi = 500)
    """
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    idx = data.obs[qc_attr].values if qc_attr is not None else np.ones(data.shape[0], dtype = bool)

    if isinstance(features, str):
        features = [features]
    if len(features) > 8:
        logger.warning("With with_control == False, only up to 8 features are allowed!")
        return None

    exprs = []
    feats = []

    size = idx.sum()
    for feature in features:
        fid = data.var_names.get_loc(feature)
        exprs.append(data.get_matrix("arcsinh.transformed")[idx, fid].toarray()[:, 0])
        feats.append(np.repeat(feature, size))

    df = pd.DataFrame({"expression": np.concatenate(exprs), "feature": np.concatenate(feats)})

    g = sns.FacetGrid(df, row="feature", hue="feature", aspect=8, height=1.0)
    g.map(sns.kdeplot, "expression", clip_on=False, shade=True, alpha=1, lw=1.5)
    g.map(sns.kdeplot, "expression", clip_on=False, color="k", lw=1)
    g.map(plt.axhline, y=0, lw=1, clip_on=False)

    def _set_label(value, color, label):
        ax = plt.gca()
        ax.text(0, 0.2, label, color="k", ha="right", va="center", transform=ax.transAxes)

    g.map(_set_label, "expression")

    g.fig.subplots_adjust(hspace=-overlap)

    g.set_titles("")
    g.set_xlabels("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    g.fig.set_dpi(dpi)
    g.fig.set_figwidth(panel_size[0])
    g.fig.set_figheight(panel_size[1])

    sns.reset_orig()

    return g.fig if return_fig else None



def doublet_plot(
    scores: List[float],
    codes: List[int],
    code_names: Optional[List[str]] = ["singlet", "singlet2", "doublet"],
    panel_size: Optional[Tuple[float, float]] = (4, 3),
    left: Optional[float] = 0.2,
    bottom: Optional[float] = 0.2,
    wspace: Optional[float] = 0.2,
    hspace: Optional[float] = 0.2,
    return_fig: Optional[bool] = False,
    dpi: Optional[float] = 300.0,
    **kwargs,
) -> Union[plt.Figure, None]:
    """Generate KDE for doublet predictions

    Parameters
    ----------
    scores: `List[float]`,
        Doublet scores.
    codes: `List[int]`
        Doublet type codes.
    code_names: `List[str]`, optional (default: `["singlet", "singlet2", "doublet"]`)
        Doublet type names.
    panel_size: `tuple`, optional (default: `(4, 3)`)
        The plot size (width, height) in inches.
    left: `float`, optional (default: `0.15`)
        This parameter sets the figure's left margin as a fraction of panel's width (left * panel_size[0]).
    bottom: `float`, optional (default: `0.15`)
        This parameter sets the figure's bottom margin as a fraction of panel's height (bottom * panel_size[1]).
    wspace: `float`, optional (default: `0.3`)
        This parameter sets the width between panels and also the figure's right margin as a fraction of panel's width (wspace * panel_size[0]).
    hspace: `float`, optional (defualt: `0.15`)
        This parameter sets the height between panels and also the figure's top margin as a fraction of panel's height (hspace * panel_size[1]).
    return_fig: ``bool``, optional, default: ``False``
        Return a ``Figure`` object if ``True``; return ``None`` otherwise.
    dpi: ``float``, optional, default: ``300.0``
        The resolution in dots per inch.

    Returns
    -------

    `Figure` object
        A ``matplotlib.figure.Figure`` object containing the dot plot if ``return_fig == True``

    Examples
    --------
    """
    from scipy.stats import gaussian_kde

    # constants
    alpha = 0.7
    eps = 1e-6

    fig, axes = _get_subplot_layouts(nrows = 2, ncols = 2, panel_size=panel_size, dpi=dpi, left=left, bottom=bottom, wspace=wspace, hspace=hspace, sharex=False, sharey=False)

    # Panel 1, KDE of scores
    ax = axes[0, 0]
    x = np.linspace(scores.min() - 0.05, scores.max() + 0.05, 200) # generate x coordinates
    kde = gaussian_kde(scores)
    y = kde(x)
    ax.plot(x, y, '-', c='k')
    ax.set_ylim(bottom = 0.0)
    ax.set_xlabel('doublet score')
    ax.set_ylabel('density')

    # Panel 2, KDE of log-transformed scores
    ax = axes[0, 1]
    scores = np.log(scores)
    x = np.linspace(scores.min() - 0.5, scores.max() + 0.5, 200) # generate x coordinates
    kde = gaussian_kde(scores)
    y = kde(x)
    ax.plot(x, y, '-', c='k')
    ax.set_ylim(bottom = 0.0)
    ax.set_xlabel('log doublet score')
    ax.set_ylabel('density')

    # Panel 3, KDE of each cluster
    ax = axes[1, 0]
    # Partition by clusters
    idx1 = codes == 0
    idx2 = codes == 1
    idx3 = codes == 2
    # Get scale factor for each cluster
    scale1 = idx1.sum() * 1.0 / codes.size
    scale2 = idx2.sum() * 1.0 / codes.size
    scale3 = idx3.sum() * 1.0 / codes.size
    # kde for singlets1
    kde1 = gaussian_kde(scores[idx1])
    y1 = kde1(x) * scale1
    idxs1 = y1 > eps
    ax.plot(x[idxs1], y1[idxs1], '-', c = 'orange', label = code_names[0])
    # kde for singlets2
    kde2 = gaussian_kde(scores[idx2])
    y2 = kde2(x) * scale2
    idxs2 = y2 > eps
    ax.plot(x[idxs2], y2[idxs2], '-', c = 'green', label = code_names[1])
    # kde for doublets
    if idx3.sum() > 0:
        scores3 = scores[idx3]
        if np.std(scores3) < eps: # if too small, add jitters
            scores3 += np.random.normal(scale = 0.01, size = scores3.size)
        kde3 = gaussian_kde(scores3)
        y3 = kde3(x) * scale3
        idxs3 = y3 > eps
        ax.plot(x[idxs3], y3[idxs3], '-', c = 'red', label = code_names[2])
    # Set ylim bottom to 0 and show legend
    ax.set_ylim(bottom = 0)
    ax.set_xlabel('log doublet score')
    ax.set_ylabel('density')
    ax.legend(loc = 'best')

    # Panel 4, overlay the overall KDE with cluster-specific KDEs
    ax = axes[1, 1]
    ax.plot(x, y, '-', c='k', label = 'overall')
    ax.plot(x[idxs1], y1[idxs1], '-', c = 'orange', alpha = alpha, label = 'singlet1')
    ax.plot(x[idxs2], y2[idxs2], '-', c = 'green', alpha = alpha, label = 'singlet2')
    if idx3.sum() > 0:
        ax.plot(x[idxs3], y3[idxs3], '-', c = 'red', alpha = alpha, label = 'doublet')
    ax.set_ylim(bottom = 0)
    ax.set_xlabel('log doublet score')
    ax.set_ylabel('density')
    ax.legend(loc = 'best')

    return fig if return_fig else None
