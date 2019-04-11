import time
import numpy as np
import igraph
from scipy.sparse import issparse



def construct_graph(W, directed = False, adjust_weights = True):
	start_time = time.time()

	assert issparse(W)

	s, t = W.nonzero()
	w = W.data
	
	if not directed:
		idx = s < t
		s = s[idx]
		t = t[idx]
		w = w[idx]

	if adjust_weights:
		w = ((w / np.median(w)) * 100.0 + 0.5).astype(int) / 100.0 # round to 2 decimal points
		idx = w > 0.0
		if idx.sum() < w.size:
			s = s[idx]
			t = t[idx]
			w = w[idx]

	G = igraph.Graph(directed = directed)
	G.add_vertices(W.shape[0])
	G.add_edges(zip(s, t))
	G.es['weight'] = w

	end_time = time.time()
	print("Graph is constructed. Time spent = {:.2f}s.".format(end_time - start_time))

	return G
