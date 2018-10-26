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

from . import row_attrs, excluded, load_10x_h5_file



def read_10x_h5_file(input_h5, genome = None, return_a_dict = False, demux_ngene = None):
	"""Load 10x-format matrices from the h5 file into a series of h5ad objects
	
	Parameters
	----------

	input_h5 : `str`
		The matricies in h5 format.
	genome : `str`, optional (default: None)
		A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
	return_a_dict : `boolean`, optional (default: False)
		If input file contains multiple genome groups, if concatenate them into one h5ad object or return a dictionary of genome-h5ad pairs. If this option is on, return a dict.
	demux_ngene : `int`, optional (default: None)
		Minimum number of genes to keep a barcode for demultiplexing.

	Returns
	-------
	
	`anndata` object or a dictionary of `anndata` objects
		An `anndata` object or a dictionary of `anndata` objects containing the count matrices.

	Examples
	--------
	>>> tools.read_10x_h5_file('example_10x.h5')
	"""	

	gdmap = load_10x_h5_file(input_h5) # gdmap , genome-data map
	if genome is not None: # remove genomes not in the list
		remove_set = set(gdmap) - set(genome.split(','))
		for gname in remove_set:
			gdmap.pop(gname)

	results = {}
	for genome, data in gdmap.items():
		# obs_dict
		barcodes = data["barcodes"].astype(str)
		if np.vectorize(lambda x: x.endswith("-1"))(barcodes).sum() == barcodes.size:
			barcodes = np.array([x[:-2] for x in barcodes]) # remove the trailing '-1'
		obs_dict = {"obs_names" : barcodes}
		obs_dict["Channel"] = ['-'.join(x.split('-')[:-1]) for x in barcodes]
		for attr, value in data.items():
			if (attr not in row_attrs) and (attr not in excluded):
				obs_dict[attr] = value.astype(str)
		# var_dict
		var_dict = {"var_names" : (data["gene_names"] if "gene_names" in data else data["antibody_names"]).astype(str)}
		if "genes" in data:
			var_dict["gene_ids"] = data["genes"].astype(str)
		# construct h5ad object
		results[genome] = anndata.AnnData(X = data["matrix"].transpose(), obs = obs_dict, var = var_dict)
		results[genome].uns["genome"] = genome

	if len(results) == 1:
		results = results[next(iter(results))]
	elif not return_a_dict:
		Xs = [] # a list of csr matrices
		vn_vec = [] # var_names vec
		gi_vec = [] # gene_ids vec
		genomes = []
		for data in results.values():
			Xs.append(data.X)
			vn_vec.append(data.var_names.values)				
			if "gene_ids" in data.var:
				gi_vec.append(data.var["gene_ids"].values)
			genomes.append(data.uns["genome"])
		var_dict = {"var_names" : np.concatenate(vn_vec)}
		if len(gi_vec) > 0:
			var_dict["gene_ids"] = np.concatenate(gi_vec)
		results = anndata.AnnData(X = hstack(Xs, format = 'csr'), obs = obs_dict, var = var_dict)
		results.uns["genome"] = ",".join(genomes)

	# for demultiplexing purpose, select only barcodes with min gene >= demux_ngene
	if demux_ngene is not None:
		assert isinstance(results, anndata.base.AnnData)
		results.obs['n_genes'] = results.X.getnnz(axis = 1)
		results.obs['n_counts'] = results.X.sum(axis = 1).A1
		results._inplace_subset_obs(results.obs['n_genes'].values >= demux_ngene)
		results.var['robust'] = True

	return results



def read_antibody_csv(input_csv):
	"""Load an ADT matrix from the csv file

	Parameters
	----------

	input_csv : `str`
		The CSV file containing ADT counts.

	Returns
	-------
	
	`anndata` object 
		An `anndata` object containing the ADT count matrix.
	
	Examples
	--------
	>>> tools.read_antibody_csv('example_ADT.csv')
	"""

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



def read_input(input_file, genome = None, return_a_dict = False, demux_ngene = None, mode = 'r+'):
	"""Load data into memory.

	This function is used to load input data into memory. Inputs can be in .h5, .h5ad, or .csv format.
	
	Parameters
	----------

	input_file : `str`
		Input file name.
	genome : `str`, .h5, optional (default: None)
		A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
	return_a_dict : `boolean`, .h5, optional (default: False)
		If input file contains multiple genome groups, if concatenate them into one h5ad object or return a dictionary of genome-h5ad pairs. If this option is on, return a dict.
	demux_ngene : `int`, .h5, optional (default: None)
		Minimum number of genes to keep a barcode for demultiplexing.
	mode : `str`, .h5ad, optional (default: `r+`)
		If input is in h5ad format, the backed mode for loading the data. mode could be 'a', 'r', 'r+'. 'a' refers to load all into memory.
	
	Returns
	-------
	`anndata` object or a dictionary of `anndata` objects
		An `anndata` object or a dictionary of `anndata` objects containing the count matrices.

	Examples
	--------
	>>> adata = tools.read_input('example_10x.h5', genome = 'mm10')
	>>> adata = tools.read_input('example.h5ad', mode = 'r+')
	>>> adata = tools.read_input('example_ADT.csv')
	"""

	if input_file.endswith('.h5'):
		data = read_10x_h5_file(input_file, genome, return_a_dict, demux_ngene)
	elif input_file.endswith('.h5ad'):
		data = anndata.read_h5ad(input_file, backed = (False if mode == 'a' else mode))
	elif input_file.endswith('.csv'):
		data = read_antibody_csv(input_file)
	else:
		print("Unrecognized file type!")
		assert False

	return data



transfer_gene_name = [(358, 'ENSG00000268991', 'FAM231C.2'), (921, 'ENSG00000278139', 'AL358075.4'), (2207, 'ENSG00000232995', 'RGS5.2'), (5847, 'ENSG00000282827', 'AC134772.2'), (5938, 'ENSG00000271858', 'CYB561D2.2'), (6087, 'ENSG00000241572', 'PRICKLE2-AS1.2'), (7213, 'ENSG00000249428', 'CFAP99.2'), (9596, 'ENSG00000280987', 'MATR3.2'), (9605, 'ENSG00000279686', 'AC142391.1'), (10277, 'ENSG00000282913', 'BLOC1S5.2'), (10867, 'ENSG00000124593', 'AL365205.1'), (11619, 'ENSG00000268592', 'RAET1E-AS1.2'), (13877, 'ENSG00000231963', 'AL662864.1'), (16117, 'ENSG00000225655', 'BX255923.1'), (16938, 'ENSG00000282955', 'RABL6.2'), (17241, 'ENSG00000265264', 'TIMM10B.2'), (18626, 'ENSG00000282682', 'C11orf71.2'), (18984, 'ENSG00000282883', 'AKR1C3.2'), (19226, 'ENSG00000150076', 'CCDC7.2'), (19346, 'ENSG00000264404', 'BX547991.1'), (21184, 'ENSG00000282031', 'TMBIM4.2'), (21230, 'ENSG00000257815', 'LINC01481.2'), (22033, 'ENSG00000228741', 'SPATA13.2'), (22037, 'ENSG00000281899', 'AL359736.3'), (22654, 'ENSG00000274827', 'LINC01297.2'), (23662, 'ENSG00000273259', 'AL049839.2'), (24019, 'ENSG00000211974', 'AC245369.1'), (26919, 'ENSG00000279257', 'C17orf100.2'), (26962, 'ENSG00000187838', 'PLSCR3'), (27137, 'ENSG00000255104', 'AC005324.4'), (27884, 'ENSG00000263715', 'LINC02210-CRHR1'), (28407, 'ENSG00000281844', 'FBF1.2'), (30440, 'ENSG00000283027', 'CAPS.2'), (32648, 'ENSG00000235271', 'LINC01422.2')]

def update_var_names(data):
	data.uns['genome'].split(',')
	prefix = re.compile('^(' + '|'.join(data.uns['genome'].split(',')) + ')_+')
	if prefix.match(data.var_names[0]):
		data.var['gene_ids'] = [prefix.sub('', x) for x in data.var['gene_ids']]
		data.var_names = pd.Index([prefix.sub('', x) for x in data.var_names])

	gsyms = data.var_names.values
	
	if data.uns['genome'] == 'GRCh38':
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
