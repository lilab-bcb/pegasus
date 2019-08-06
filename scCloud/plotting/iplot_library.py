import plotly.offline as py
from natsort import natsorted

from .plot_utils import get_palettes

import numpy as np
import pandas as pd
import scipy.sparse


def scatter(df, output_file):
    labels = natsorted(df["Annotation"].unique())
    palettes = get_palettes(len(labels))
    data = []
    for i, (label, color) in enumerate(zip(labels, palettes)):
        tdf = df[df["Annotation"] == label]
        trace = dict(
            name=label,
            x=tdf.iloc[:, 1],
            y=tdf.iloc[:, 2],
            type="scattergl",
            mode="markers",
            marker=dict(size=5, color=color, line=dict(width=0)),
        )
        data.append(trace)
    layout = dict(
        title=output_file,
        xaxis=dict(title=df.columns[1]),
        yaxis=dict(title=df.columns[2]),
        margin=dict(l=0, r=0, b=0, t=0),
    )
    fig = dict(data=data, layout=layout)
    py.plot(fig, filename=output_file)


def scatter_real(df, output_file, log10=False):
    if log10:
        df["Annotation"] = np.log10(df["Annotation"])

    trace = dict(
        x=df.iloc[:, 1],
        y=df.iloc[:, 2],
        type="scattergl",
        mode="markers",
        marker=dict(
            size=5,
            color=df["Annotation"],
            colorscale="Jet",
            colorbar=dict(title="Density"),
            line=dict(width=0),
        ),
    )
    data = [trace]
    layout = dict(
        title=output_file,
        xaxis=dict(title=df.columns[1]),
        yaxis=dict(title=df.columns[2]),
        margin=dict(l=0, r=0, b=0, t=0),
    )
    fig = dict(data=data, layout=layout)
    py.plot(fig, filename=output_file)


def scatter3d(df, output_file):
    labels = natsorted(df["Annotation"].unique())
    palettes = get_palettes(len(labels))
    data = []
    for i, (label, color) in enumerate(zip(labels, palettes)):
        tdf = df[df["Annotation"] == label]
        trace = dict(
            name=label,
            x=tdf.iloc[:, 1],
            y=tdf.iloc[:, 2],
            z=tdf.iloc[:, 3],
            type="scatter3d",
            mode="markers",
            marker=dict(size=5, color=color, line=dict(width=0)),
        )
        data.append(trace)
    layout = dict(
        title=output_file,
        scene=dict(
            xaxis=dict(title=df.columns[1]),
            yaxis=dict(title=df.columns[2]),
            zaxis=dict(title=df.columns[3]),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
    )
    fig = dict(data=data, layout=layout)
    py.plot(fig, filename=output_file)


def scatter3d_real(df, output_file, log10=False):
    if log10:
        df["Annotation"] = np.log10(df["Annotation"])

    trace = dict(
        x=df.iloc[:, 1],
        y=df.iloc[:, 2],
        z=df.iloc[:, 3],
        type="scatter3d",
        mode="markers",
        marker=dict(
            size=5,
            color=df["Annotation"],
            colorscale="Jet",
            colorbar=dict(title="Density"),
            line=dict(width=0),
        ),
    )
    data = [trace]
    layout = dict(
        title=output_file,
        scene=dict(
            xaxis=dict(title=df.columns[1]),
            yaxis=dict(title=df.columns[2]),
            zaxis=dict(title=df.columns[3]),
        ),
        margin=dict(l=0, r=0, b=0, t=0),
    )
    fig = dict(data=data, layout=layout)
    py.plot(fig, filename=output_file)
