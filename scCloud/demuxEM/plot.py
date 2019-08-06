import numpy as np
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from natsort import natsorted


def plot_adt_hist(adt, attr, out_file, alpha=0.5, dpi=500, figsize=None):
    idx_signal = np.isin(adt.obs[attr], "signal")
    signal = adt.obs.loc[idx_signal, "counts"]
    background = adt.obs.loc[~idx_signal, "counts"]
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
    data, out_file, plot_attr="n_counts", cat_attr="demux_type", dpi=500, figsize=None
):
    bins = np.logspace(
        np.log10(min(data.obs[plot_attr])), np.log10(max(data.obs[plot_attr])), 101
    )
    cat_vec = data.obs[cat_attr]
    ax = plt.gca()
    if cat_attr == "demux_type":
        ax.hist(
            data.obs.loc[np.isin(cat_vec, "singlet"), plot_attr],
            bins,
            alpha=0.5,
            label="singlet",
        )
        ax.hist(
            data.obs.loc[np.isin(cat_vec, "doublet"), plot_attr],
            bins,
            alpha=0.5,
            label="doublet",
        )
        ax.hist(
            data.obs.loc[np.isin(cat_vec, "unknown"), plot_attr],
            bins,
            alpha=0.5,
            label="unknown",
        )
    ax.legend(loc="upper right")
    ax.set_xscale("log")
    ax.set_xlabel("Number of RNA UMIs (log10 scale)")
    ax.set_ylabel("Number of cellular barcodes")
    if figsize is not None:
        plt.gcf().set_size_inches(*figsize)
    plt.savefig(out_file, dpi=dpi)
    plt.close()


def plot_bar(heights, tick_labels, xlabel, ylabel, out_file, dpi=500, figsize=None):
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


def plot_dataframe_bar(df, ylabel, out_file, dpi=500, figsize=None):
    if df.shape[1] == 1:
        df.plot.bar(legend=False)
    else:
        df.plot.bar()
    ax = plt.gca()
    ax.set_ylabel(ylabel)
    if figsize is not None:
        plt.gcf().set_size_inches(*figsize)
    plt.savefig(out_file, dpi=dpi)
    plt.close()


# attrs is a dict with name: attr format; if this is a gene violin, attrs == {gene: gene_name}
def plot_violin(
    data,
    attrs,
    out_file,
    xlabel=None,
    ylabel=None,
    title=None,
    dpi=500,
    figsize=None,
    linewidth=None,
    log=False,
    inner="box",
):
    df = None

    if "gene" in attrs:
        df = pd.DataFrame(
            data[:, attrs["gene"]].X.toarray(),
            index=data.obs_names,
            columns=[attrs["gene"]],
        )
        df["assignment"] = data.obs["demux_type"].astype(str)
        idx_singlet = np.isin(data.obs["demux_type"], "singlet")
        singlets = data.obs.loc[idx_singlet, "assignment"].astype(str)
        df.loc[idx_singlet, "assignment"] = singlets
        categories = natsorted(singlets.unique())
        categories.extend(["doublet", "unknown"])
        df["assignment"] = pd.Categorical(df["assignment"], categories=categories)
        xlabel = "assignment"
        ylabel = attrs["gene"]
    else:
        dfs = []
        if isinstance(data, anndata.base.AnnData):
            for name, attr in attrs.items():
                dfs.append(pd.DataFrame({xlabel: name, ylabel: data.obs[attr].values}))
        else:
            for arr, name in zip(data, attrs):
                dfs.append(pd.DataFrame({xlabel: name, ylabel: arr}))
        df = pd.concat(dfs)

    if log:
        df[ylabel] = np.log10(df[ylabel])
        sns.violinplot(
            x=xlabel, y=ylabel, data=df, linewidth=linewidth, cut=0, inner=inner
        )
        y_max = int(np.ceil(df[ylabel].max()))
        loc = list(range(y_max + 1))
        labels = [r"$10^{}$".format(x) for x in loc]
        plt.yticks(loc, labels)
    else:
        sns.violinplot(
            x=xlabel, y=ylabel, data=df, linewidth=linewidth, cut=0, inner=inner
        )

    ax = plt.gca()
    ax.grid(False)
    if "gene" in attrs:
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


def plot_heatmap(vec1, vec2, out_file, dpi=500, xlabel="", ylabel=""):
    df = pd.crosstab(vec1, vec2)
    df.columns.name = df.index.name = ""

    ax = plt.gca()
    ax.xaxis.tick_top()
    ax = sns.heatmap(df, annot=True, fmt="d", cmap="inferno", ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.tight_layout()
    plt.savefig(out_file, dpi=500)
    plt.close()


# type = count or gene
def plot_human_vs_mouse(
    data, out_file, plot_type, alpha=0.5, dpi=500, log=False, figsize=None
):
    type2xlabel = {
        "count": "Number of mouse RNA UMIs",
        "gene": "Number of expressed mouse genes",
    }
    type2ylabel = {
        "count": "Number of human RNA UMIs",
        "gene": "Number of expressed human genes",
    }
    type2xattr = {"count": "mm10_n_counts", "gene": "mm10_n_genes"}
    type2yattr = {"count": "grch38_n_counts", "gene": "grch38_n_genes"}

    ax = plt.gca()

    labels = data.obs["demux_type"].astype(str)
    labels[
        np.isin(
            data.obs["assignment"],
            ["Sample1MmF", "Sample2MmF", "Sample3MmM", "Sample4MmM"],
        )
    ] = "singlet_mouse"
    # labels[np.isin(data.obs['assignment'], ['Sample5HuF','Sample6HuF','Sample7HuM','Sample8HuM'])] = 'singlet_human'
    labels[
        np.isin(data.obs["assignment"], ["Sample5HuF", "Sample6HuF", "Sample7HuM"])
    ] = "singlet_human"
    labels[np.isin(data.obs["assignment"], "Sample8HuM")] = "singlet_donor8"

    idx_dbl = np.isin(labels, "doublet")
    mouse_set = {"Sample1MmF", "Sample2MmF", "Sample3MmM", "Sample4MmM"}
    # human_set = {'Sample5HuF','Sample6HuF','Sample7HuM','Sample8HuM'}
    human_set = {"Sample5HuF", "Sample6HuF", "Sample7HuM"}

    def doublet_type(assign_str):
        ids = set(assign_str.split(","))
        # return 'doublet_human' if ids.issubset(human_set) else ('doublet_mouse' if ids.issubset(mouse_set) else 'doublet_mix')
        return (
            "doublet_human"
            if ids.issubset(human_set)
            else (
                "doublet_mouse"
                if ids.issubset(mouse_set)
                else ("doublet_donor8" if "Sample8HuM" in ids else "doublet_mix")
            )
        )

    labels[idx_dbl] = np.vectorize(doublet_type)(
        data.obs.loc[idx_dbl, "assignment"].values
    )

    # labels = pd.Categorical(labels, categories = ['singlet_human', 'singlet_mouse', 'doublet_human', 'doublet_mouse', 'doublet_mix', 'unknown'])
    # colors = ['red', 'blue', 'purple', 'green', 'fuchsia', 'grey']
    labels = pd.Categorical(
        labels,
        categories=[
            "singlet_human",
            "singlet_mouse",
            "singlet_donor8",
            "doublet_human",
            "doublet_mouse",
            "doublet_mix",
            "doublet_donor8",
            "unknown",
        ],
    )
    colors = [
        "red",
        "blue",
        "greenyellow",
        "purple",
        "green",
        "fuchsia",
        "orange",
        "grey",
    ]

    for k, cat in enumerate(labels.categories):
        idx = np.isin(labels, cat)
        ax.scatter(
            data.obs.loc[idx, type2xattr[plot_type]],
            data.obs.loc[idx, type2yattr[plot_type]],
            c=colors[k],
            marker=".",
            alpha=alpha,
            edgecolors="none",
            label=cat,
            rasterized=True,
        )

    ax.grid(False)
    ax.legend()
    ax.set_xlabel(type2xlabel[plot_type])
    ax.set_ylabel(type2ylabel[plot_type])

    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")

    if figsize is not None:
        plt.gcf().set_size_inches(*figsize)
    plt.savefig(out_file, dpi=dpi)
    plt.close()
