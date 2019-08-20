import time
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import igraph
from collections import deque 

from typing import List



def calc_pseudotime(data: "AnnData", roots: List[str]) -> None:
    """
    TODO: documentation.
    """
    start = time.time()

    if not isinstance(roots, list):
        roots = [roots]

    if "X_diffmap" not in data.obsm.keys():
        raise ValueError("Please run diffmap first!")

    data.uns["roots"] = roots
    mask = np.isin(data.obs_names, data.uns["roots"])
    distances = np.mean(
        euclidean_distances(data.obsm["X_diffmap"][mask, :], data.obsm["X_diffmap"]),
        axis=0,
    )
    dmin = distances.min()
    dmax = distances.max()
    data.obs["pseudotime"] = (distances - dmin) / (dmax - dmin)

    end = time.time()
    print(
        "calc_pseudotime finished. Time spent = {:.2f}s".format(end - start)
    )



def construct_knn_graph(indices, distances):
    G = igraph.Graph(directed=False)
    G.add_vertices(indices.shape[0])
    edges = []
    w = []
    for i in range(indices.shape[0]):
        for j in range(indices.shape[1]):
            edges.append((i, indices[i][j]))
            w.append(distances[i][j])
    G.add_edges(edges)
    G.es["weight"] = w
    return G

def bfs_on_mst(G, root_id):
    mst = G.spanning_tree(weights = "weight")
    myiter = mst.bfsiter(root_id, advanced = True)
    n = G.vcount()
    parents = np.full(n, -1, dtype = int)
    for value in myiter:
        if value[2] is not None:
            parents[value[0].index] = value[2].index
    return parents

def infer_path(data, cluster, clust_id, path_name, k: int = 10):
    assert 'roots' in data.uns and len(data.uns['roots']) == 1
    root_id = int(np.isin(data.obs_names, data.uns['roots'][0]).nonzero()[0][0])
    indices = data.uns['diffmap_knn_indices']
    distances = data.uns['diffmap_knn_distances']
    G = construct_knn_graph(indices, distances)    
    parents = bfs_on_mst(G, root_id)
    inpath = np.zeros(data.shape[0], dtype = bool)
    idx = np.isin(data.obs[cluster], clust_id)
    inpath[idx] = True

    qsize = idx.sum()
    queue = deque(idx.nonzero()[0])
    while qsize > 0:
        vid = queue.popleft()
        qsize -= 1
        if parents[vid] >= 0 and not inpath[parents[vid]]:
            inpath[parents[vid]] = True
            queue.append(parents[vid])
            qsize += 1

    for vid in np.nonzero(inpath & ~idx)[0]:
        inpath[indices[vid, 0:k]] = True

    data.obs[path_name] = inpath.astype(str)





