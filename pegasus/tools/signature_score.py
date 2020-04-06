import numpy as np
import pandas as pd

from typing import Dict, List
from anndata import AnnData

import logging

logger = logging.getLogger("pegasus")
from pegasus.utils import decorators as pg_deco


@pg_deco.TimeLogger()
def calc_signature_score(data: AnnData, signatures: Dict[str, List[str]], n_bins: int = 50):
	if "mean" not in data.var:
		data.var["mean"] = data.X.mean(axis = 0).A1

	if data.uns.get("sig_n_bins", 0) != n_bins:
		data.uns["sig_n_bins"] = n_bins
		mean_vec = data.var["mean"]
		bins = pd.qcut(mean_vec, n_bins)
		bins.cat.categories = bins.cat.categories.astype(str)
		data.var["bins"] = bins
		# calculate background expectations
		sig_background = np.zeros((data.shape[0], n_bins))
		for code in range(n_bins):
			idx = (bins.cat.codes == code).values
			base = mean_vec[idx].mean()
			sig_background[:, code] = data.X[:, idx].mean(axis = 1).A1 - base
		data.obsm["sig_background"] = sig_background

	for key, gene_list in signatures.items():
		genes = pd.Index(gene_list)
		idx = data.var_names.isin(genes)
		if idx.sum() < genes.size:
			omitted = ~genes.isin(data.var_names)
			logger.warning("For signature {}, genes {} are not in the data and thus omitted!".format(key, str(list(genes[omitted]))[1:-1]))
		if key in data.obs:
			logger.warning("Signature key {} exists in data.obs, the existing content will be overwritten!".format(key))
		data.obs[key] = (data.X[:, idx].mean(axis = 1).A1 - data.var.loc[idx, "mean"].mean()) - data.obsm["sig_background"][:, data.var["bins"].cat.codes[idx]].mean(axis = 1)
