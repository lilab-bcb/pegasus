import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from natsort import natsorted

from typing import Tuple, List
from pegasusio import UnimodalData
from plotting import scatter


def plot_hto_hist(hashing_data: UnimodalData, attr: str, out_file: str, alpha: float = 0.5, dpi: int = 500, figsize: Tuple[float, float] = None) -> None:
    idx_signal = np.isin(hashing_data.obs[attr], "signal")
    signal = hashing_data.obs.loc[idx_signal, "counts"]
    background = hashing_data.obs.loc[~idx_signal, "counts"]
    bins = np.logspace(0, np.log10(max(signal.max(), background.max())), 501)
    plt.hist(background, bins, alpha=alpha, label="background", log=True)
    plt.hist(signal, bins, alpha=alpha, label="signal", log=True)
    plt.legend(loc="upper right")
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_xlabel("Number of hashtag UMIs (log10 scale)")
    ax.set_ylabel("Number of cellular barcodes (log10 scale)")
    if figsize is not None:
        plt.gcf().set_size_inches(*figsize)
    plt.savefig(out_file, dpi=dpi)
    plt.close()


def plot_rna_hist(
    rna_data: UnimodalData, out_file: str, plot_attr: str = "n_counts", cat_attr: str = "demux_type", dpi: int = 500, figsize: Tuple[float, float] = None
) -> None:
    fig, axes = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 9], 'hspace': 0.6}, figsize = figsize, dpi = dpi)

    # Percentage of RNA barcodes having HTO tags
    nhto = rna_data.obs_names.sum()
    total = rna_data.shape[0]

    ax = axes[0]
    p1 = nhto * 100.0 / total
    p2 = (total - nhto) * 100.0 / total
    labels = ['Has HTO', 'No HTO']
    ax.barh(y = 0, width = p1, color = 'red', height = 0.5, label = labels[0])
    ax.barh(y = 0, width = p2, left = p1, color = 'gray', height = 0.5, label = labels[1])
    ax.set_yticks([])
    ax.set_ylabel('RNA    \nbarcodes', rotation = 0, ha = 'right', va = 'center')
    ax.set_xlim(left = 0, right = 100)
    ax.set_xlabel('Percentage')
    ax.legend(ncol=2, bbox_to_anchor=(1, 1), loc='lower right', fontsize = 'small')

    # RNA histogram
    bins = np.logspace(
        np.log10(min(rna_data.obs[plot_attr])), np.log10(max(rna_data.obs[plot_attr])), 101
    )
    cat_vec = rna_data.obs[cat_attr]
    ax = axes[1]
    ax.hist(
        rna_data.obs.loc[np.isin(cat_vec, "singlet"), plot_attr],
        bins,
        alpha=0.5,
        label="singlet",
    )
    ax.hist(
        rna_data.obs.loc[np.isin(cat_vec, "doublet"), plot_attr],
        bins,
        alpha=0.5,
        label="doublet",
    )
    ax.hist(
        rna_data.obs.loc[np.isin(cat_vec, "unknown"), plot_attr],
        bins,
        alpha=0.5,
        label="unknown",
    )
    ax.legend(loc="upper right")
    ax.set_xscale("log")
    ax.set_xlabel("Number of RNA UMIs (log10 scale)")
    ax.set_ylabel("Number of cellular barcodes")

    fig.savefig(out_file)


def plot_hto_background(hashing_data: UnimodalData, xlabel: str, ylabel: str, out_file: str, dpi: int = 500, figsize: Tuple[float, float] = None) -> None:
    
    heights = hashing_data.uns["background_probs"]
    tick_labels = hashing_data.var_names,
    plt.bar(
        x=np.linspace(0.5, heights.size - 0.5, heights.size),
        height=heights,
        tick_label=tick_labels,
    )
    ax = plt.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if figsize is not None:
        plt.gcf().set_size_inches(*figsize)
    rotation = 90 if max([len(x) for x in tick_labels]) > 6 else 0
    plt.tick_params(axis="x", labelsize=7, labelrotation=rotation)
    plt.tight_layout()
    plt.savefig(out_file, dpi=dpi)
    plt.close()

def plot_gene_violin(
    data: UnimodalData,
    gene_name: str,
    out_file: str,
    title: str = None,
    dpi: int = 500,
    figsize: Tuple[float, float] = None,
    linewidth: float = None,
    inner: str = "box",
) -> None:
    df = pd.DataFrame(
        data[:, gene_name].X.toarray(),
        index=data.obs_names,
        columns=[gene_name],
    )
    df["assignment"] = data.obs["demux_type"].astype(str)
    idx_singlet = np.isin(data.obs["demux_type"], "singlet")
    singlets = data.obs.loc[idx_singlet, "assignment"].astype(str)
    df.loc[idx_singlet, "assignment"] = singlets
    categories = natsorted(singlets.unique())
    categories.extend(["doublet", "unknown"])
    df["assignment"] = pd.Categorical(df["assignment"], categories=categories)
    xlabel = "assignment"
    ylabel = gene_name

    sns.violinplot(
        x=xlabel, y=ylabel, data=df, linewidth=linewidth, cut=0, inner=inner
    )

    ax = plt.gca()
    ax.grid(False)
    ax.set_ylabel("log(TP100K+1)")
    if title is not None:
        ax.set_title(title)

    if figsize is not None:
        plt.gcf().set_size_inches(*figsize)

    rotation = 90 if max([len(x) for x in df[xlabel].unique()]) > 6 else 0
    plt.tick_params(axis="x", labelsize=7, labelrotation=rotation)
    plt.tight_layout()
    plt.savefig(out_file, dpi=dpi)
    plt.close()

def plot_hto_umap(
    data: UnimodalData,
    out_file: str,
    dpi: int = 500,
) -> None: 
    fig = scatter(data, basis='umap', dpi=dpi, return_fig=True)
    fig.savefig(out_file)