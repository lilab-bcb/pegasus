import time
import numpy as np 
import pandas as pd
import anndata
import tables
import xlsxwriter

from scipy.sparse import issparse, csr_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.utils.extmath import randomized_svd



obs_not = set(["data", "indices", "indptr", "shape", "gene_names", "genes", "barcodes"])

def read_10x_h5_file(input_h5, genome):
	with tables.open_file(input_h5) as h5_in:
		inpmat = {}
		for node in h5_in.walk_nodes("/" + genome, "Array"):
			inpmat[node.name] = node.read()

	X = csr_matrix((inpmat["data"], inpmat["indices"], inpmat["indptr"]), shape = (inpmat["shape"][1], inpmat["shape"][0]))
	
	obs_dict = {"obs_names" : inpmat["barcodes"].astype(str)}
	for key, value in inpmat.items():
		if key not in obs_not:
			obs_dict[key] = value.astype(str)

	data = anndata.AnnData(X = X, obs = obs_dict, var = {"var_names" : inpmat["gene_names"].astype(str), "gene_ids": inpmat["genes"].astype(str)})

	return data

def read_input(input_name, genome = 'GRCh38', is_10x = True):
	data = read_10x_h5_file(input_name + '_10x.h5', genome) if is_10x else anndata.read_h5ad(input_name + '.h5ad')
	data.obs['Channel'] = ['-'.join(x.split('-')[:-2]) for x in data.obs_names]
	return data



transfer_gene_name = [(358, 'ENSG00000268991', 'FAM231C.2'), (921, 'ENSG00000278139', 'AL358075.4'), (2207, 'ENSG00000232995', 'RGS5.2'), (5847, 'ENSG00000282827', 'AC134772.2'), (5938, 'ENSG00000271858', 'CYB561D2.2'), (6087, 'ENSG00000241572', 'PRICKLE2-AS1.2'), (7213, 'ENSG00000249428', 'CFAP99.2'), (9596, 'ENSG00000280987', 'MATR3.2'), (9605, 'ENSG00000279686', 'AC142391.1'), (10277, 'ENSG00000282913', 'BLOC1S5.2'), (10867, 'ENSG00000124593', 'AL365205.1'), (11619, 'ENSG00000268592', 'RAET1E-AS1.2'), (13877, 'ENSG00000231963', 'AL662864.1'), (16117, 'ENSG00000225655', 'BX255923.1'), (16938, 'ENSG00000282955', 'RABL6.2'), (17241, 'ENSG00000265264', 'TIMM10B.2'), (18626, 'ENSG00000282682', 'C11orf71.2'), (18984, 'ENSG00000282883', 'AKR1C3.2'), (19226, 'ENSG00000150076', 'CCDC7.2'), (19346, 'ENSG00000264404', 'BX547991.1'), (21184, 'ENSG00000282031', 'TMBIM4.2'), (21230, 'ENSG00000257815', 'LINC01481.2'), (22033, 'ENSG00000228741', 'SPATA13.2'), (22037, 'ENSG00000281899', 'AL359736.3'), (22654, 'ENSG00000274827', 'LINC01297.2'), (23662, 'ENSG00000273259', 'AL049839.2'), (24019, 'ENSG00000211974', 'AC245369.1'), (26919, 'ENSG00000279257', 'C17orf100.2'), (26962, 'ENSG00000187838', 'PLSCR3'), (27137, 'ENSG00000255104', 'AC005324.4'), (27884, 'ENSG00000263715', 'LINC02210-CRHR1'), (28407, 'ENSG00000281844', 'FBF1.2'), (30440, 'ENSG00000283027', 'CAPS.2'), (32648, 'ENSG00000235271', 'LINC01422.2')]

def update_var_names(data, genome):
	if data.var_names[0].startswith(genome + "_"):
		n = len(genome) + 1
		data.var['gene_ids'] = [x[n:] for x in data.var['gene_ids']]
		data.var_names = pd.Index([x[n:] for x in data.var_names])

	gsyms = data.var_names.values
	
	if genome == "GRCh38":
		for pos, gid, gsym in transfer_gene_name:
			assert data.var.iloc[pos, 0] == gid
			gsyms[pos] = gsym
	else:	
		dup_ids = Counter()
		for i in range(gsyms.size):
			idn = dup_ids[gsyms[i]]
			dup_ids[gsyms[i]] += 1
			if idn > 0:
				gsyms[i] = gsyms[i] + ".{}".format(idn)
	
	data.var_names = pd.Index(gsyms)



def filter_data(data, mito_prefix = 'MT-', filt_xlsx = None, min_genes = 500, max_genes = 6000, percent_mito = 0.1, percent_cells = 0.0005):
	data.obs['n_genes'] = data.X.getnnz(axis = 1)
	data.obs['n_counts'] = data.X.sum(axis = 1).A1
	mito_genes = [name for name in data.var_names if name.startswith(mito_prefix)]
	data.obs['percent_mito'] = data[:, mito_genes].X.sum(axis=1).A1 / data.obs['n_counts'].values

	# Filter cells	
	if filt_xlsx is not None:
		writer = pd.ExcelWriter(filt_xlsx, engine='xlsxwriter')
		tot_c = data.obs['Channel'].value_counts()
	
	obs_index = np.logical_and.reduce((data.obs['n_genes'] >= min_genes, 
									   data.obs['n_genes'] < max_genes,
									   data.obs['percent_mito'] < percent_mito))
	data._inplace_subset_obs(obs_index)

	if filt_xlsx is not None:
		kept_c = data.obs['Channel'].value_counts().reindex(tot_c.index, fill_value = 0)
		df = pd.DataFrame({'Kept' : kept_c, 'Filt' : tot_c - kept_c, 'Total' : tot_c})
		df = df[['Kept', 'Filt', 'Total']]
		df.sort_values('Kept', inplace = True)
		df.to_excel(writer, sheet_name = "Cell filtration stats")

	# Filter genes
	data.var['n_cells'] = data.X.getnnz(axis = 0)
	data.var['percent_cells'] = data.var['n_cells'] / data.shape[0]
	data.var['robust'] = data.var['percent_cells'] >= percent_cells

	if filt_xlsx is not None:
		idx = data.var['robust'] == False
		df = pd.DataFrame({'n_cells': data.var.loc[idx, 'n_cells'], 'percent_cells': data.var.loc[idx, 'percent_cells']})
		df.sort_values('n_cells', ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Gene filtration stats")
		writer.save()
		print("Filtration results are written.")

	var_index = (data.var['n_cells'] > 0).values
	data._inplace_subset_var(var_index)
	print("After filteration, {nc} cells and {ng} genes are kept. Among {ng} genes, {nrb} genes are robust.".format(nc = data.shape[0], ng = data.shape[1], nrb = data.var['robust'].sum()))


def log_norm(data, norm_count):
	""" Normalization and then take log """
	assert issparse(data.X)
	mat = data.X[:, data.var['robust'].values]
	scale = norm_count / mat.sum(axis = 1).A1
	data.X.data *= np.repeat(scale, np.diff(data.X.indptr))
	data.X = data.X.log1p()

def run_pca(data, standardize = True, max_value = 10, nPC = 50, random_state = 0):
	start = time.time()
	if issparse(data.X):
		data.X = data.X.toarray()

	if standardize:
		scaler = StandardScaler(copy = False)
		scaler.fit_transform(data.X)

	if max_value is not None:
		data.X[data.X > max_value] = max_value

	pca = PCA(n_components = nPC, random_state = random_state)
	X_pca = pca.fit_transform(data.X)	
	data.obsm['X_pca'] = X_pca
	data.varm['PCs'] = pca.components_.T
	data.uns['pca'] = {}
	data.uns['pca']['variance'] = pca.explained_variance_
	data.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_
	end = time.time()
	print("PCA is done. Time spent = {:.2f}s.".format(end - start))

def run_rpca(data, scale = False, max_value = 10.0, nPC = 50, random_state = 0):
	""" smooth outliers, then no center/scale data """
	start = time.time()

	# Smooth out outliers
	means, variances = mean_variance_axis(data.X, axis = 0)
	stds = np.sqrt(variances * (data.X.shape[0] / (data.X.shape[0] - 1))) # make it unbiased
	assert (stds == 0.0).sum() == 0

	data_new = (data.X.data - means[data.X.indices]) / stds[data.X.indices]
	outliers = data_new > max_value
	data.X.data[outliers] = max_value * stds[data.X.indices[outliers]] + means[data.X.indices[outliers]]

	if scale:
		data.X.data /= stds[data.X.indices]

	U, S, VT = randomized_svd(data.X, n_components = nPC, random_state = random_state)
	data.obsm['X_rpca'] = U * S

	end = time.time()
	print("RPCA is done. Time spent = {:.2f}s.".format(end - start))
