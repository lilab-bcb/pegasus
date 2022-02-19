import time
import numpy as np
from scipy.sparse import issparse, csr_matrix

import logging
logger = logging.getLogger(__name__)
try:
    import igraph
except ImportError:
    logger.error("Need python-igraph!  Try 'pip install python-igraph'.")

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


@timer(logger=logger)
def construct_graph(
    W: csr_matrix, directed: bool = False, adjust_weights: bool = True
) -> "igraph":

    assert issparse(W)

    s, t = W.nonzero()
    w = W.data

    if not directed:
        idx = s < t
        s = s[idx]
        t = t[idx]
        w = w[idx]

    if adjust_weights:
        w = ((w / np.median(w)) * 100.0 + 0.5).astype(
            int
        ) / 100.0  # round to 2 decimal points
        idx = w > 0.0
        if idx.sum() < w.size:
            s = s[idx]
            t = t[idx]
            w = w[idx]

    G = igraph.Graph(directed=directed)
    G.add_vertices(W.shape[0])
    G.add_edges(zip(s, t))
    G.es["weight"] = w

    return G
