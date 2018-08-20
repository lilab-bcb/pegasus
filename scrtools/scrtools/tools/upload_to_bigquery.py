import re
import numpy as np
import pandas as pd
from google.cloud import bigquery
import threading
import time


def make_it_unique(names):
	pattern = re.compile('\W')
	seen = set()
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

def get_metadata_data_frames(data):
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

	df2 = pd.DataFrame({'name' : ['PCA_X', 'PCA_Y', 'TSNE_X', 'TSNE_Y', 'DIFFMAP_X', 'DIFFMAP_Y', 'DIFFMAP_Z'], 
						'min_value' : [data.obsm[pca_str][:,0].min(), data.obsm[pca_str][:, 1].min(), data.obsm['X_tsne'][:, 0].min(), data.obsm['X_tsne'][:, 1].min(), data.obsm['X_diffmap_pca'][:, 0].min(), data.obsm['X_diffmap_pca'][:, 1].min(), data.obsm['X_diffmap_pca'][:, 2].min()],
						'max_value' : [data.obsm[pca_str][:,0].max(), data.obsm[pca_str][:, 1].max(), data.obsm['X_tsne'][:, 0].max(), data.obsm['X_tsne'][:, 1].max(), data.obsm['X_diffmap_pca'][:, 0].max(), data.obsm['X_diffmap_pca'][:, 1].max(), data.obsm['X_diffmap_pca'][:, 2].max()]})

	return df, df2

def upload_table(df, client, dataset, table):
	start = time.time()

	table_ref = client.dataset(dataset).table(table)

	try:
		client.delete_table(table_ref)
	except:
		None

	job = client.load_table_from_dataframe(df, table_ref)
	job.result()

	end = time.time()
	print("Table {}.{} is uploaded. Time spent = {:.2f}s.".format(dataset, table, end - start))

	return job.state == 'DONE'

def upload_expression_matrix(thread_no, dfs_genes, n_threads, batch_size, X, gene_names, bigquery_names, metadata_df, client, dataset, table_prefix):
	n_genes = X.shape[1]
	for i, start_pos in enumerate(range(batch_size * thread_no, n_genes, batch_size * n_threads)):
		table_id = thread_no + i * n_threads
		end_pos = min(start_pos + batch_size, n_genes)
		dfs_genes[table_id] = pd.DataFrame({'gene' : gene_names[start_pos : end_pos], 'bigquery_name' : bigquery_names[start_pos : end_pos], 'table_number' : table_id})
		df_expr = pd.DataFrame(data = X[:, start_pos : end_pos].toarray(), columns = bigquery_names[start_pos : end_pos])
		df_expr[metadata_df.columns] = metadata_df
		upload_table(df_expr, client, dataset, table_prefix + '_expr{}'.format(table_id))

def upload_data_to_bigquery(data, client, dataset, table_prefix, batch_size = 1000, n_threads = -1):
	start = time.time()

	metadata_df, range_df = get_metadata_data_frames(data)
	X = data.X.tocsc()
	gene_names = data.var_names.values
	bigquery_names = make_it_unique(gene_names)

	n_tables = (X.shape[1] // batch_size) + (X.shape[1] % batch_size > 0)
	if n_threads < 0:
		n_threads = n_tables

	threads = [None] * n_threads
	dfs_genes = [None] * n_tables
	for i in range(n_threads):
		t = threading.Thread(target=upload_expression_matrix, args=(i, dfs_genes, n_threads, batch_size, X, gene_names, bigquery_names, metadata_df, client, dataset, table_prefix))
		threads[i] = t
		t.start()

	for i in range(n_threads):
		threads[i].join()

	df_genes = pd.concat(dfs_genes)
	upload_table(df_genes, client, dataset, table_prefix + '_gene_names')
	upload_table(range_df, client, dataset, table_prefix + '_ranges')

	end = time.time()
	print("Data were uploaded to BigQuery, total time = {:.2f}s.".format(end - start))

