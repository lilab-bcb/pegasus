import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from typing import List, Tuple


def _transform_basis(basis: str) -> str:
    if basis == "tsne" or basis == "citeseq_tsne":
        return "tSNE"
    elif basis == "fitsne":
        return "FItSNE"
    elif basis == "umap":
        return "UMAP"
    elif basis == "diffmap":
        return "DC"
    elif basis == "pca" or basis == "rpca":
        return "PC"
    elif basis == "diffmap_pca":
        return "DPC"
    elif basis == "fle":
        return "FLE"
    elif basis == "net_tsne":
        return "Net-tSNE"
    elif basis == "net_fitsne":
        return "Net-FItSNE"
    elif basis == "net_umap":
        return "Net-UMAP"
    elif basis == "net_fle":
        return "Net-FLE"
    else:
        return basis

def _get_dot_size(size_arr, size_min, size_max, dot_min, dot_max):
    size_pixel = np.interp(size_arr, (size_min, size_max), (dot_min, dot_max))
    size_pixel = 5 * size_pixel
    return size_pixel

def _get_nrows_and_ncols(num_figs: int, nrows: int, ncols: int) -> Tuple[int, int]:
    if nrows is None and ncols is None:
        ncols = int(np.sqrt(num_figs - 1)) + 1
        nrows = (num_figs // ncols) + (num_figs % ncols > 0)
    elif nrows is None:
        nrows = (num_figs // ncols) + (num_figs % ncols > 0)
    elif ncols is None:
        ncols = (num_figs // nrows) + (num_figs % nrows > 0)

    return nrows, ncols


def _get_marker_size(nsamples: int) -> float:
    return min(20.0, (240000.0 if nsamples > 300000 else 120000.0) / nsamples)


def _get_subplot_layouts(
    nrows: int = 1,
    ncols: int = 1,
    panel_size: Tuple[float, float] = (6, 4),
    dpi: float = 300.0,
    left: float = 0.0,
    bottom: float = 0.0,
    wspace: float = 0.0,
    hspace: float = 0.0,
    squeeze: bool = True,
    sharex: bool = True,
    sharey: bool = True,
    frameon: bool = False,
):
    """
    Get subplot layouts. Default is nrows = ncols = 1
    """
    left_margin = left * panel_size[0]
    bottom_margin = bottom * panel_size[1]
    right_space = wspace * panel_size[0]
    top_space = hspace * panel_size[1]

    figsize = (
        left_margin + panel_size[0] * (1.0 + wspace) * ncols,
        bottom_margin + panel_size[1] * (1.0 + hspace) * nrows,
    )
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        squeeze=squeeze,
        sharex=sharex,
        sharey=sharey,
        frameon=frameon,
        dpi = dpi,
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


def _get_legend_ncol(label_size: int, max_ncol: int = None):
    max_ncol = 100 if max_ncol is None else max_ncol
    return min(1 if label_size <= 14 else (2 if label_size <= 30 else 3), max_ncol)





# palettes are imported from scanpy, need to be replaced

pegasus_20 = [
    "#c5b0d5",
    "#ff7f0e",
    "#8c564b",
    "#ff9896",
    "#1f77b4",
    "#dbdb8d",
    "#e377c2",
    "#2ca02c",
    "#aec7e8",
    "#ffbb78",
    "#9edae5",
    "#98df8a",
    "#d62728",
    "#9467bd",
    "#c49c94",
    "#f7b6d2",
    "#bcbd22",
    "#17becf",
    "#ad494a",
    "#8c6d31",
]

zeileis_26 = [
    "#023fa5",
    "#7d87b9",
    "#bec1d4",
    "#d6bcc0",
    "#bb7784",
    "#8e063b",
    "#4a6fe3",
    "#8595e1",
    "#b5bbe3",
    "#e6afb9",
    "#e07b91",
    "#d33f6a",
    "#11c638",
    "#8dd593",
    "#c6dec7",
    "#ead3c6",
    "#f0b98d",
    "#ef9708",
    "#0fcfc0",
    "#9cded6",
    "#d5eae7",
    "#f3e1eb",
    "#f6c4e1",
    "#f79cd4",
    "#7f7f7f",
    "#c7c7c7",
    "#1CE6FF",
    "#336600",  # these last ones were added,
]

godsnot_64 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00",
    "#1CE6FF",
    "#FF34FF",
    "#FF4A46",
    "#008941",
    "#006FA6",
    "#A30059",
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
    "#5A0007",
    "#809693",
    "#FEFFE6",
    "#1B4400",
    "#4FC601",
    "#3B5DFF",
    "#4A3B53",
    "#FF2F80",
    "#61615A",
    "#BA0900",
    "#6B7900",
    "#00C2A0",
    "#FFAA92",
    "#FF90C9",
    "#B903AA",
    "#D16100",
    "#DDEFFF",
    "#000035",
    "#7B4F4B",
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
    "#6A3A4C",
]


def _get_palettes(n_labels: int, with_background: bool = False, show_background: bool = False):
    if with_background:
        n_labels -= 1

    if n_labels <= 20:
        palettes = pegasus_20
    elif n_labels <= 26:
        palettes = zeileis_26
    else:
        assert n_labels <= 64
        palettes = godsnot_64

    if with_background:
        palettes = np.array(
            ["gainsboro" if show_background else "white"] + palettes[:n_labels]
        )
    else:
        palettes = np.array(palettes[:n_labels])

    return palettes

def _plot_cluster_labels_in_heatmap(ax, ids, orientation):
    pos_acc = 0
    ticks = []
    labels = []
    max_label_len = 0
    for label, value in ids.value_counts().iteritems():
        labels.append(label)
        if len(label) > max_label_len:
            max_label_len = len(label)
        ticks.append(pos_acc + value / 2)
        pos_acc += value

    ax.grid(False)

    if orientation == 'left':
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, rotation=0, x=-0.05)
        ax.tick_params(axis='y', left=False)
    else:
        color_box = ax.get_position()
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, rotation=90, y=color_box.y1+max_label_len/100)
        ax.xaxis.set_ticks_position('top')
        ax.tick_params(axis='x', top=False)

    # Remove borders
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


Restriction = namedtuple("Restriction", ["negation", "values"])


class RestrictionParser:
    def __init__(self, restrictions: List[str]):
        self.restrs = {}

        if restrictions is None:
            return None

        for restr_str in restrictions:
            attr, value_str = restr_str.split(":")
            negation = False
            if value_str[0] == "~":
                negation = True
                value_str = value_str[1:]
            self.restrs[attr] = Restriction(negation=negation, values=value_str.split(","))


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
