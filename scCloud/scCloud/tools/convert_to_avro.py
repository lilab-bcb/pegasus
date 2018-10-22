import numpy as np
import fastavro 
import time


def generate_schema(data, name):
	schema = {
		'doc' : 'Single cell data in avro format.',
		'name' : name,
		'namespace' : 'scCloud',
		'type' : 'record'
	}

	fields = [{'name' : 'barcode', 'type' : 'string'}]
	
	for value in data.obs.columns:
		col = {'name' : value}
		dtype = data.obs.dtypes[value]
		if dtype == np.bool:
			dtype = 'boolean'
		elif dtype == np.int32:
			dtype = 'int'
		elif dtype == np.int64:
			dtype = 'long'
		elif dtype == np.float32:
			dtype = 'float'
		elif dtype == np.float64:
			dtype = 'double'
		else:
			dtype = 'string'
		col['type'] = dtype
		fields.append(col)
	
	for gene_name in data.var_names:
		fields.append({'name' : gene_name, 'type' : 'float', 'default' : 0.0})
	
	schema['fields'] = fields
	
	return schema



def generate_records(data):
	time_a = time.time()

	pca_str = 'X_pca' if 'X_pca' in data.obsm.keys() else 'X_rpca'
	X_pca = data.obsm[pca_str][:, 0:2]
	X_tsne = data.obsm['X_tsne']
	X_diffmap = data.obsm['X_diffmap_pca']

	data.obs.index.name = 'barcode'
	df = data.obs.reset_index()

	genes = data.var_names.values
	indptr = data.X.indptr
	indices = data.X.indices
	expr = data.X.data

	time_b = time.time()
	print("Preprocess time = {:.2f}".format(time_b - time_a))
	time_a = time_b

	nsample = data.shape[0]
	for i in range(nsample):
		rec_dict = {genes[y] : expr[indptr[i] + x] for x, y in enumerate(indices[indptr[i] : indptr[i + 1]])}
		rec_dict.update(df.iloc[i].to_dict())
		rec_dict['PCA_X'] = X_pca[i, 0]
		rec_dict['PCA_Y'] = X_pca[i, 1]
		rec_dict['TSNE_X'] = X_tsne[i, 0]
		rec_dict['TSNE_Y'] = X_tsne[i, 1]
		rec_dict['DIFFMAP_X'] = X_diffmap[i, 0]
		rec_dict['DIFFMAP_Y'] = X_diffmap[i, 1]
		rec_dict['DIFFMAP_Z'] = X_diffmap[i, 2]
		yield rec_dict
		if (i + 1) % 1000 == 0:
			time_b = time.time()
			print("Processed {} records, time spent = {:.2f}.".format(i + 1, time_b - time_a))
			time_a = time_b

	print("Time spent = {:.2f}".format(time_b - time_a))


def write_data_to_avro(data, name, outfile):
	schema = generate_schema(data, name)
	parsed_schema = fastavro.parse_schema(schema)
	with open(outfile, 'wb') as out:
		fastavro.writer(out, parsed_schema, generate_records(data), codec = 'snappy', validator = True)

