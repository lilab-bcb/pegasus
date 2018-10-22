import re
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import time


seen = set()
pattern = re.compile('\W')

def make_it_unique(names):
	global seen

	results = []
	for name in names:
		new_name = pattern.sub('_', name)
		stub = new_name.lower()
		suf = ''
		num = 1
		while (stub + suf) in seen:
			num += 1
			suf = '_' + str(num)
		seen.add(stub + suf)	
		new_name = 'gene_' + new_name + suf
		results.append(new_name)

	return results



def convert_to_parquet(data, output_name, max_num_gene = 500):
	time_start = time.time()

	data.obs.index.name = 'barcode'
	df = data.obs.reset_index()
	pca_str = 'X_pca' if 'X_pca' in data.obsm.keys() else 'X_rpca'
	df['PCA_X'] = data.obsm[pca_str][:, 0]
	df['PCA_Y'] = data.obsm[pca_str][:, 1]
	df['TSNE_X'] = data.obsm['X_tsne'][:, 0]
	df['TSNE_Y'] = data.obsm['X_tsne'][:, 1]
	df['DIFFMAP_X'] = data.obsm['X_diffmap_pca'][:, 0]
	df['DIFFMAP_Y'] = data.obsm['X_diffmap_pca'][:, 1]
	df['DIFFMAP_Z'] = data.obsm['X_diffmap_pca'][:, 2]

	metadata_table = pa.Table.from_pandas(df)

	df_genes = pd.DataFrame(columns = ['gene', 'bigquery_name', 'table_number'])
	X = data.X.tocsc()

	for i, start in enumerate(range(0, data.shape[1], max_num_gene)):
		end = min(start + max_num_gene, data.shape[1])

		output_file = output_name + '.expr.{}.parquet'.format(i)
		bigquery_names = make_it_unique(data.var_names[start:end])

		df_genes.append(pd.DataFrame({'gene' : data.var_names[start:end], 'bigquery_name' : bigquery_names, 'table_number' : i}))

		df_expr = pd.DataFrame(data = X[:, start:end].toarray(), columns = bigquery_names)
		table = pa.Table.from_pandas(df_expr)
		for column in metadata_table.itercolumns():
			if column.name != '__index_level_0__':
				table = table.append_column(column)

		pq.write_table(table, output_file)
		print(output_file + ' is written!')

	output_file = output_name + '.gene_name_table.parquet'
	pq.write_table(pa.Table.from_pandas(df_genes), output_file)
	print(output_file + ' is written!')

	time_end = time.time()
	print("Total time spent = {:.2f}s.".format(time_end - time_start))
