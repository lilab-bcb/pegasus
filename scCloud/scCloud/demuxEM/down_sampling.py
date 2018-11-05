import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import estimate_background_probs, demultiplex
from ..tools import read_input


def down_sampling(data_gt, adt_gt, probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], n_threads = 1):
	f = np.vectorize(lambda x, p: np.random.binomial(int(x + 1e-4), p, size = 1)[0])

	nsample = data_gt.shape[0]
	nhto = adt_gt.X.sum()
	
	fracs = []	
	accuracy = []
	for p in probs:		
		data = data_gt.copy()
		adt = adt_gt.copy()
		
		adt.X.data = f(adt.X.data, p)
		idx = adt.X.sum(axis = 1).A1 > 0
		adt = adt[idx, ].copy()
		fracs.append(adt.X.sum() / nhto)

		estimate_background_probs(adt)
		demultiplex(data, adt, n_threads = n_threads)
		accuracy.append(sum(data.obs['assignment'].values.astype('str') == data_gt.obs['assignment'].values.astype('str')) * 100.0 / nsample)

	fracs.append(1.0)
	accuracy.append(100.0)

	return fracs, accuracy


def plot_down_sampling(rna_file, adt_file, out_file, probs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], n_threads = 1, dpi = 500, figsize = None):
	data_gt = read_input(rna_file, mode = 'a')
	adt_gt = read_input(adt_file, mode = 'a')
	fracs, accuracy = down_sampling(data_gt, adt_gt, probs = probs, n_threads = n_threads)
	plt.plot(fracs, accuracy, 'o-')
	ax = plt.gca()
	ax.set_xlabel("Fraction of HTO UMIs")
	ax.set_ylabel("Accuracy (in percentage)")
	if figsize is not None:
		plt.gcf().set_size_inches(*figsize)
	plt.savefig(out_file, dpi = dpi)
	plt.close()
