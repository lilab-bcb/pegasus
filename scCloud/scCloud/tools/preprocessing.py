import time
import re
import numpy as np 
import pandas as pd
import anndata
import tables
import xlsxwriter
from collections import Counter, ChainMap

from scipy.sparse import issparse, csr_matrix, hstack
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.utils.extmath import randomized_svd



def read_10x_h5_file(input_h5, genome, ngene = None):
	obs_not = set(["data", "indices", "indptr", "shape", "gene_names", "genes", "barcodes"])

	Xs = [] # a list of csr matrices
	gn_vec = [] # gene_names vec
	gi_vec = [] # gene_ids vec

	obs_dict = None

	with tables.open_file(input_h5) as h5_in:
		genomes = []
		if genome is None: # if no genome is provided, scan the hdf5 file, must load raw matrix here.
			for i, group in enumerate(h5_in.walk_groups()):
				if i > 0:
					genomes.append(group._v_name)
		else:
			genomes.append(genome)

		for genome in genomes:
			inpmat = {}
			for node in h5_in.walk_nodes("/" + genome, "Array"):
				inpmat[node.name] = node.read()
			
			Xs.append(csr_matrix((inpmat["data"], inpmat["indices"], inpmat["indptr"]), shape = (inpmat["shape"][1], inpmat["shape"][0])))
			gn_vec.append(inpmat["gene_names"].astype(str))
			gi_vec.append(inpmat["genes"].astype(str))

			if obs_dict is None:
				barcodes = inpmat["barcodes"].astype(str)
				if np.vectorize(lambda x: x.endswith("-1"))(barcodes).sum() == barcodes.size:
					barcodes = [x[:-2] for x in barcodes] # remove the trailing '-1'
				obs_dict = {"obs_names" : barcodes}
				obs_dict["Channel"] = ['-'.join(x.split('-')[:-1]) for x in barcodes]
				for key, value in inpmat.items():
					if key not in obs_not:
						obs_dict[key] = value.astype(str)

	data = anndata.AnnData(X = Xs[0] if len(Xs) == 1 else hstack(Xs, format = 'csr'),
						   obs = obs_dict, 
						   var = {"var_names" : gn_vec[0] if len(gn_vec) == 1 else np.concatenate(gn_vec),
						   		  "gene_ids" : gi_vec[0] if len(gi_vec) == 1 else np.concatenate(gi_vec)})

	if ngene is not None:
		n_genes_vec = data.X.getnnz(axis = 1)
		data._inplace_subset_obs(n_genes_vec >= ngene)

	return data

def read_antibody_file(input_h5at, genome):
	obs_not = set(["data", "indices", "indptr", "shape", "antibody_names", "barcodes"])

	with tables.open_file(input_h5at) as h5_in:
		inpmat = {}
		for node in h5_in.walk_nodes("/" + genome, "Array"):
			inpmat[node.name] = node.read()

	X = csr_matrix((inpmat["data"], inpmat["indices"], inpmat["indptr"]), shape = (inpmat["shape"][1], inpmat["shape"][0]))
	
	barcodes = inpmat["barcodes"].astype(str)
	obs_dict = {"obs_names" : barcodes}
	obs_dict["Channel"] = ['-'.join(x.split('-')[:-1]) for x in barcodes]
	for key, value in inpmat.items():
		if key not in obs_not:
			obs_dict[key] = value.astype(str)

	data = anndata.AnnData(X = X, obs = obs_dict, var = {"var_names" : inpmat["antibody_names"].astype(str)})

	return data

def read_antibody_csv(input_csv):
	barcodes = []
	antibody_names = []
	stacks = []
	with open(input_csv) as fin:
		barcodes = next(fin).strip().split(',')[1:]
		for line in fin:
			fields = line.strip().split(',')
			antibody_names.append(fields[0])
			stacks.append([int(x) for x in fields[1:]])
	data = anndata.AnnData(X = csr_matrix(np.stack(stacks, axis = 1)), obs = {"obs_names" : barcodes}, var = {"var_names" : antibody_names})

	return data

def read_input(input_file, genome = None, mode = 'r+', ngene = None):
	"""Load either 10x-formatted raw count matrix or h5ad-formatted processed expression matrix into memory.

	This function is used to load input data into memory.
	
	Parameters
	----------

	input_file : `str`
		Input file name.
	genome : `str`, optional (default: None)
		The genome used to produce raw count matrices. If genome == None, we will load count matrices from all possible genomes and merge them into one big matrix.
	mode : `str`, optional (default: `r+`)
		If input is h5ad format, the backed mode for loading the data. mode could be 'a', 'r', 'r+'. 'a' refers to load all into memory.
	ngene : `int`, optional (default: None)
		Only used for raw 10x hdf5 file. If set, only keep cells/nuclei with at least <ngene> expressed genes.
	Returns
	-------
	`anndata` object
		An `anndata` object contains the count matrix.

	Examples
	--------
	>>> adata = tools.read_input('example_10x.h5', genome = 'mm10')
	>>> adata = tools.read_input('example.h5ad', mode = 'r+')
	"""

	if input_file.endswith('.h5'):
		data = read_10x_h5_file(input_file, genome, ngene = ngene)
	elif input_file.endswith('.h5ad'):
		data = anndata.read_h5ad(input_file, backed = (False if mode == 'a' else mode))
	elif input_file.endswith('.h5at'):
		data = read_antibody_file(input_file, genome)
	elif input_file.endswith('.csv'):
		data = read_antibody_csv(input_file)
	else:
		print("Unrecognized file type!")
		assert False

	return data



transfer_gene_name = [(358, 'ENSG00000268991', 'FAM231C.2'), (921, 'ENSG00000278139', 'AL358075.4'), (2207, 'ENSG00000232995', 'RGS5.2'), (5847, 'ENSG00000282827', 'AC134772.2'), (5938, 'ENSG00000271858', 'CYB561D2.2'), (6087, 'ENSG00000241572', 'PRICKLE2-AS1.2'), (7213, 'ENSG00000249428', 'CFAP99.2'), (9596, 'ENSG00000280987', 'MATR3.2'), (9605, 'ENSG00000279686', 'AC142391.1'), (10277, 'ENSG00000282913', 'BLOC1S5.2'), (10867, 'ENSG00000124593', 'AL365205.1'), (11619, 'ENSG00000268592', 'RAET1E-AS1.2'), (13877, 'ENSG00000231963', 'AL662864.1'), (16117, 'ENSG00000225655', 'BX255923.1'), (16938, 'ENSG00000282955', 'RABL6.2'), (17241, 'ENSG00000265264', 'TIMM10B.2'), (18626, 'ENSG00000282682', 'C11orf71.2'), (18984, 'ENSG00000282883', 'AKR1C3.2'), (19226, 'ENSG00000150076', 'CCDC7.2'), (19346, 'ENSG00000264404', 'BX547991.1'), (21184, 'ENSG00000282031', 'TMBIM4.2'), (21230, 'ENSG00000257815', 'LINC01481.2'), (22033, 'ENSG00000228741', 'SPATA13.2'), (22037, 'ENSG00000281899', 'AL359736.3'), (22654, 'ENSG00000274827', 'LINC01297.2'), (23662, 'ENSG00000273259', 'AL049839.2'), (24019, 'ENSG00000211974', 'AC245369.1'), (26919, 'ENSG00000279257', 'C17orf100.2'), (26962, 'ENSG00000187838', 'PLSCR3'), (27137, 'ENSG00000255104', 'AC005324.4'), (27884, 'ENSG00000263715', 'LINC02210-CRHR1'), (28407, 'ENSG00000281844', 'FBF1.2'), (30440, 'ENSG00000283027', 'CAPS.2'), (32648, 'ENSG00000235271', 'LINC01422.2')]

def update_var_names(data, genome):
	prefix = re.compile('^' + genome + '_+')
	if prefix.match(data.var_names[0]):
		data.var['gene_ids'] = [prefix.sub('', x) for x in data.var['gene_ids']]
		data.var_names = pd.Index([prefix.sub('', x) for x in data.var_names])

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
	data.obs['percent_mito'] = data[:, mito_genes].X.sum(axis=1).A1 / np.maximum(data.obs['n_counts'].values, 1.0)

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


def filter_cells_cite_seq(data, max_cells):
	data.obs['n_counts'] = data.X.sum(axis = 1).A1
	obs_index = np.zeros(data.shape[0], dtype = bool)
	obs_index[np.argsort(data.obs['n_counts'].values)[::-1][:max_cells]] = True
	data._inplace_subset_obs(obs_index)
	data.var['robust'] = True
	print("After filteration, {nc} cells are kept, with the minimum nUMI = {numi}.".format(nc = max_cells, numi = data.obs['n_counts'].min()))

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

def parse_subset_selections(subset_selections):	
	subsets_dict = {}
	for subset_str in subset_selections:
		print(subset_str)
		attr, value_str = subset_str.split(':')
		if attr in subsets_dict:
			subsets_dict[attr].extend(value_str.split(','))
		else:
			subsets_dict[attr] = value_str.split(',')
	return subsets_dict

def get_anndata_for_subclustering(data, subset_selections):	
	obs_index = np.full_like(data.obs[data.obs.columns.values[0]], True)
	subsets_dict = parse_subset_selections(subset_selections)
	for key, value in subsets_dict.items():
		print(key, 'corresponds to', subsets_dict[key])
		obs_index = obs_index & np.isin(data.obs[key], value)	
	data = data[obs_index, :]
	obs_dict = {"obs_names" : data.obs_names.values}
	for attr in data.obs.columns:
		if attr != "pseudotime":
			if attr.find("_labels") < 0:
				obs_dict[attr] = data.obs[attr].values
			else:
				obs_dict['parent_' + attr] = data.obs[attr].values

	var_dict = {"var_names" : data.var_names.values, "gene_ids": data.var["gene_ids"].values, "robust": data.var["robust"]}

	newdata = anndata.AnnData(X = data.X, obs = obs_dict, var = var_dict)
	if "Channels" in data.uns:
		newdata.uns["Channels"] = data.uns["Channels"]
	if "Groups" in data.uns:
		newdata.uns["Groups"] = data.uns["Groups"]
	if "means" in data.varm.keys():
		newdata.varm["means"] = data.varm["means"]
	if "stds" in data.varm.keys():
		newdata.varm["stds"] = data.varm["stds"]

	print("{0} cells are selected.".format(newdata.shape[0]))
	
	return newdata
