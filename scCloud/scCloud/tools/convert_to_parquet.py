import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import time
from . import read_input



def convert_to_parquet(data, output_name, nthreads):
	data.obs.index.name = 'barcode'
	df = data.obs.reset_index()
	if 'X_pca' in data.obsm.keys():
		df['PCA_X'] = data.obsm['X_pca'][:, 0]
		df['PCA_Y'] = data.obsm['X_pca'][:, 1]
	if 'X_rpca' in data.obsm.keys():
		df['RPCA_X'] = data.obsm['X_rpca'][:, 0]
		df['RPCA_Y'] = data.obsm['X_rpca'][:, 1]
	if 'X_tsne' in data.obsm.keys():
		df['TSNE_X'] = data.obsm['X_tsne'][:, 0]
		df['TSNE_Y'] = data.obsm['X_tsne'][:, 1]
	if 'X_fitsne' in data.obsm.keys():
		df['FITSNE_X'] = data.obsm['X_fitsne'][:, 0]
		df['FITSNE_Y'] = data.obsm['X_fitsne'][:, 1]
	if 'X_umap' in data.obsm.keys():
		df['UMAP_X'] = data.obsm['X_umap'][:, 0]
		df['UMAP_Y'] = data.obsm['X_umap'][:, 1]
	if 'X_fle' in data.obsm.keys():
		df['FLE_X'] = data.obsm['X_fle'][:, 0]
		df['FLE_Y'] = data.obsm['X_fle'][:, 1]
	if 'X_diffmap_pca' in data.obsm.keys():		
		df['DIFFMAP_X'] = data.obsm['X_diffmap_pca'][:, 0]
		df['DIFFMAP_Y'] = data.obsm['X_diffmap_pca'][:, 1]
		df['DIFFMAP_Z'] = data.obsm['X_diffmap_pca'][:, 2]
	metadata_table = pa.Table.from_pandas(df, nthreads = nthreads)


	df_expr = pd.DataFrame(data = data.X.toarray(), columns = data.var_names)
	parquet_table = pa.Table.from_pandas(df_expr, nthreads = nthreads)

	for i in range(metadata_table.num_columns - 1, 0, -1):
		column = metadata_table[i]
		if column.name != '__index_level_0__':
			parquet_table = parquet_table.add_column(0, column)

	output_file = output_name + '.parquet'
	pq.write_table(parquet_table, output_file)
	print(output_file + ' is written!')



def run_conversion(input_h5ad_file, output_name, nthreads):
	start = time.time()
	data = read_input(input_h5ad_file, mode = 'a')
	end = time.time()
	print("Time spent for loading the expression matrix is {:.2f}s.".format(end - start))
	
	start = time.time()
	convert_to_parquet(data, output_name, nthreads)
	end = time.time()
	print("Time spent on generating the PARQUET file is {:.2f}s.".format(end - start))
