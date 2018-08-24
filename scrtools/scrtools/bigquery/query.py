import numpy as np
import pandas as pd
import uuid
import time
import threading
import tempfile
from google.cloud import bigquery, storage, exceptions

class QueryEngine:

	MAX_PAGE_SIZE = 100000 # 100K
	MIN_PAGE_SIZE = 50000 # 50K
	NUM_BINS_2D = 1024 # default 2D discretization granule
	NUM_BINS_3D = 512 # default 3D discretization granule
	NUM_THREADS = 10 # default number of threads to query BigQuery
	MAX_SAMPLE_SIZE_NODISC = 1000000 # maximum sample size to not to discretize 

	QUERY_TEMPLATES = {
		'gene_names' : 'SELECT * FROM `{dataset}.gene_names`',
		'ranges' : 'SELECT * FROM `{dataset}.ranges`',

		'2D_real' : 'SELECT ROUND({basis}_X, 5) AS {basis}_X, ROUND({basis}_Y, 5) AS {basis}_Y, ROUND({attr}, 2) AS {attr_name} FROM `{dataset}.expression_{table_id}`',
		'2D_cat' : 'SELECT ROUND({basis}_X, 5) AS {basis}_X, ROUND({basis}_Y, 5) AS {basis}_Y, {attr} AS {attr_name} FROM `{dataset}.expression_{table_id}`',
		'3D_real' : 'SELECT ROUND({basis}_X, 5) AS {basis}_X, ROUND({basis}_Y, 5) AS {basis}_Y, ROUND({basis}_Z, 5) AS {basis}_Z, ROUND({attr}, 2) AS {attr_name} FROM `{dataset}.expression_{table_id}`',
		'3D_cat' : 'SELECT ROUND({basis}_X, 5) AS {basis}_X, ROUND({basis}_Y, 5) AS {basis}_Y, ROUND({basis}_Z, 5) AS {basis}_Z, {attr} AS {attr_name} FROM `{dataset}.expression_{table_id}`',

		'2D_count' : 'SELECT COUNT(DISTINCT CAST(FLOOR(({basis}_X - {x_min}) / {x_denom} * {bin_size}) AS INT64) * {bin_size} + CAST(FLOOR(({basis}_Y - {y_min}) / {y_denom} * {bin_size}) AS INT64)) FROM `{dataset}.expression_{table_id}`',
		'3D_count' : 'SELECT COUNT(DISTINCT CAST(FLOOR(({basis}_X - {x_min}) / {x_denom} * {bin_size}) AS INT64) * {bin_size} * {bin_size} + CAST(FLOOR(({basis}_Y - {y_min}) / {y_denom} * {bin_size}) AS INT64) * {bin_size} + CAST(FLOOR(({basis}_Z - {z_min}) / {z_denom} * {bin_size}) AS INT64)) FROM `{dataset}.expression_{table_id}`',
		
		'2D_real_disc' : ('SELECT ROUND(AVG({basis}_X), 5) AS {basis}_X, '
								 'ROUND(AVG({basis}_Y), 5) AS {basis}_Y, '
								 'ROUND(AVG({attr}), 2) AS {attr_name}, '
								 'CAST(FLOOR(({basis}_X - {x_min}) / {x_denom} * {bin_size}) AS INT64) * {bin_size} + CAST(FLOOR(({basis}_Y - {y_min}) / {y_denom} * {bin_size}) AS INT64) AS MYID '
						  'FROM `{dataset}.expression_{table_id}` '
						  'GROUP BY MYID'),
		'2D_cat_disc' : ('SELECT ROUND(AVG({basis}_X), 5) AS {basis}_X, '
								'ROUND(AVG({basis}_Y), 5) AS {basis}_Y, '
								'APPROX_TOP_COUNT({attr}, 1)[OFFSET(0)].value AS {attr_name}, '
								'CAST(FLOOR(({basis}_X - {x_min}) / {x_denom} * {bin_size}) AS INT64) * {bin_size} + CAST(FLOOR(({basis}_Y - {y_min}) / {y_denom} * {bin_size}) AS INT64) AS MYID '
						 'FROM `{dataset}.expression_{table_id}` '
						 'GROUP BY MYID'),

		'3D_real_disc' : ('SELECT ROUND(AVG({basis}_X), 5) AS {basis}_X, '
								 'ROUND(AVG({basis}_Y), 5) AS {basis}_Y, '
								 'ROUND(AVG({basis}_Z), 5) AS {basis}_Z, '
								 'ROUND(AVG({attr}), 2) AS {attr_name}, '
								 'CAST(FLOOR(({basis}_X - {x_min}) / {x_denom} * {bin_size}) AS INT64) * {bin_size} * {bin_size} + CAST(FLOOR(({basis}_Y - {y_min}) / {y_denom} * {bin_size}) AS INT64) * {bin_size} + CAST(FLOOR(({basis}_Z - {z_min}) / {z_denom} * {bin_size}) AS INT64) AS MYID '
						  'FROM `{dataset}.expression_{table_id}` '
						  'GROUP BY MYID'),
		'3D_cat_disc' : ('SELECT ROUND(AVG({basis}_X), 5) AS {basis}_X, '
								'ROUND(AVG({basis}_Y), 5) AS {basis}_Y, '
								'ROUND(AVG({basis}_Z), 5) AS {basis}_Z, '
								'APPROX_TOP_COUNT({attr}, 1)[OFFSET(0)].value AS {attr_name}, '
								'CAST(FLOOR(({basis}_X - {x_min}) / {x_denom} * {bin_size}) AS INT64) * {bin_size} * {bin_size} + CAST(FLOOR(({basis}_Y - {y_min}) / {y_denom} * {bin_size}) AS INT64) * {bin_size} + CAST(FLOOR(({basis}_Z - {z_min}) / {z_denom} * {bin_size}) AS INT64) AS MYID '
						 'FROM `{dataset}.expression_{table_id}` '
						 'GROUP BY MYID')
	}



	def __init__(self, project_id, dataset, is_parallel = True, location = 'US'):
		self.project_id = project_id
		self.dataset = dataset
		self.is_parallel = is_parallel
		self.location = location
		
		self.client = bigquery.Client(project = project_id)

		if not is_parallel:
			self.bucket_id = 'temp-vis-' + uuid.uuid4().hex
			storage_client = storage.Client(project = project_id)
			try:
				self.bucket = storage_client.get_bucket(self.bucket_id)
			except exceptions.NotFound:
				self.bucket = storage_client.create_bucket(self.bucket_id)

		# Query gene names
		start_time = time.time()
		query = QueryEngine.QUERY_TEMPLATES['gene_names'].format(dataset = self.dataset)
		query_job = self.client.query(query, location = self.location)
		df_tmp = query_job.result().to_dataframe()
		self.df_gene_names = df_tmp[['bigquery_name', 'table_number']]
		self.df_gene_names.index = df_tmp['gene']
		end_time = time.time()
		print("Time = {:.2f}s.".format(end_time - start_time))

		# Query ranges
		start_time = time.time()
		query = QueryEngine.QUERY_TEMPLATES['ranges'].format(dataset = self.dataset)
		query_job = self.client.query(query, location = self.location)
		df_tmp = query_job.result().to_dataframe()
		self.df_ranges = df_tmp[['min_value', 'max_value']]
		self.df_ranges.index = df_tmp['name']
		end_time = time.time()
		print("Time = {:.2f}s.".format(end_time - start_time))

		# Parse table schema
		table = self.client.get_table(self.client.dataset(self.dataset).table('expression_0'))
		self.nsample = table.num_rows
		self.attr2type = {}
		for schema_field in table.schema:
			attr = schema_field.name
			if not attr.startswith('gene_') and not attr.endswith('_X') and not attr.endswith('_Y') and not attr.endswith('_Z') and attr != '__index_level_0__':
				self.attr2type[schema_field.name] = schema_field.field_type

		# discretization dict from (basis, bin_size) to num_rows
		self.disc_dict = {}



	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		if not self.is_parallel and exc_type is None: # exit successfully
			self.bucket.delete(force = True)



	def calc_num_rows(self, basis, bin_size):
		num_rows = self.disc_dict.get((basis, bin_size), None)
		if num_rows is None:
			query_string = ''
			if basis != 'DIFFMAP':
				query_string = QueryEngine.QUERY_TEMPLATES['2D_count'].format(basis = basis, bin_size = bin_size, dataset = self.dataset, table_id = 0,
					x_min = self.df_ranges.at[basis + '_X', 'min_value'], x_denom = self.df_ranges.at[basis + '_X', 'max_value'] - self.df_ranges.at[basis + '_X', 'min_value'],
					y_min = self.df_ranges.at[basis + '_Y', 'min_value'], y_denom = self.df_ranges.at[basis + '_Y', 'max_value'] - self.df_ranges.at[basis + '_Y', 'min_value'])
			else:
				query_string = QueryEngine.QUERY_TEMPLATES['3D_count'].format(basis = basis, bin_size = bin_size, dataset = self.dataset, table_id = 0,
					x_min = self.df_ranges.at[basis + '_X', 'min_value'], x_denom = self.df_ranges.at[basis + '_X', 'max_value'] - self.df_ranges.at[basis + '_X', 'min_value'],
					y_min = self.df_ranges.at[basis + '_Y', 'min_value'], y_denom = self.df_ranges.at[basis + '_Y', 'max_value'] - self.df_ranges.at[basis + '_Y', 'min_value'],
					z_min = self.df_ranges.at[basis + '_Z', 'min_value'], z_denom = self.df_ranges.at[basis + '_Z', 'max_value'] - self.df_ranges.at[basis + '_Z', 'min_value'])
			print(query_string)
			start_time = time.time()
			query_job = self.client.query(query_string, location = self.location)
			num_rows = list(query_job.result())[0].values()[0]
			end_time = time.time()
			print("Time calc num rows ({basis}, {bin_size}) = {time:.2f}.".format(basis = basis, bin_size = bin_size, time = end_time - start_time))	

			self.disc_dict[basis, bin_size] = num_rows

		return num_rows



	def query_storage(self, query):
		
		# Note, transferring to google bucket could be slow occasionally
		start_time = time.time()
		query_job = self.client.query(query, location = self.location)
		query_job.result()
		end_time = time.time()
		print("Time querying = {:.2f}s.".format(end_time - start_time))

		start_time = end_time
		table_file = "{}.csv.gz".format(query_job.destination.table_id)
		config = bigquery.ExtractJobConfig()
		config.compression = 'GZIP' # GZIP is faster in transferring table and downloading table locally
		extract_job = self.client.extract_table(query_job.destination, 'gs://{bucket}/{table_file}'.format(bucket = self.bucket.id, table_file = table_file), job_config = config)
		extract_job.result()
		end_time = time.time()
		print("Time transfering = {:.2f}s.".format(end_time - start_time))

		start_time = end_time
		blob = self.bucket.get_blob(table_file)
		fp = tempfile.TemporaryFile()
		blob.download_to_file(fp)
		fp.seek(0)
		df = pd.read_csv(fp, compression = 'gzip')
		end_time = time.time()
		print("Time downloading = {:.2f}s.".format(end_time - start_time))

		return df



	def submit_query(self, thread_no, query_string, page_size, df_arr):
		start_time = time.time()

		query = query_string + ' LIMIT {count} OFFSET {skip_rows}'.format(count = page_size, skip_rows = page_size * thread_no)
		query_job = self.client.query(query, location = self.location)
		df_arr[thread_no] = query_job.result().to_dataframe()

		end_time = time.time()
		# print("Thread = {}, Time = {:.2f}.".format(thread_no, end_time - start_time))


	def query_parallel(self, query_string, num_rows):
		start_time = time.time()
		
		num_threads = QueryEngine.NUM_THREADS
		page_size = max((num_rows // num_threads) + (num_rows % num_threads > 0), QueryEngine.MIN_PAGE_SIZE)
		if page_size > QueryEngine.MAX_PAGE_SIZE:
			print("Warning: page size = {} > {}, query time might be long!".format(page_size, QueryEngine.MAX_PAGE_SIZE))
		num_threads = (num_rows // page_size) + (num_rows % page_size > 0)
		print("num_threads = {}, page_size = {}".format(num_threads, page_size))

		threads = [None] * num_threads
		df_arr = [None] * num_threads
		for i in range(num_threads):
			t = threading.Thread(target = self.submit_query, args = (i, query_string, page_size, df_arr))
			threads[i] = t
			t.start()

		for i in range(num_threads):
			threads[i].join()

		end_time = time.time()
		print("Time spent = {:.2f}s".format(end_time - start_time))

		return pd.concat(df_arr)



	def query_gene(self, basis, gene_name, discretization = True, bin_size = None):
		results = None

		if gene_name not in self.df_gene_names.index:
			print("Query failed: {gene_name} does not exist!".format(gene_name = gene_name))
			return results

		bigquery_name = self.df_gene_names.at[gene_name, 'bigquery_name']
		table_id = self.df_gene_names.at[gene_name, 'table_number']

		# If sample size > MAX_SAMPLE_SIZE_NODISC, set discretization to True
		if not discretization and self.nsample > QueryEngine.MAX_SAMPLE_SIZE_NODISC:
			discretization = True

		query_string = ''
		num_rows = 0
		if discretization:
			if basis != 'DIFFMAP':
				if bin_size is None:
					bin_size = QueryEngine.NUM_BINS_2D
				query_string = QueryEngine.QUERY_TEMPLATES['2D_real_disc'].format(basis = basis, attr = bigquery_name, attr_name = gene_name, 
					x_min = self.df_ranges.at[basis + '_X', 'min_value'], x_denom = self.df_ranges.at[basis + '_X', 'max_value'] - self.df_ranges.at[basis + '_X', 'min_value'],
					y_min = self.df_ranges.at[basis + '_Y', 'min_value'], y_denom = self.df_ranges.at[basis + '_Y', 'max_value'] - self.df_ranges.at[basis + '_Y', 'min_value'],
					bin_size = bin_size, dataset = self.dataset, table_id = table_id)
			else:
				if bin_size is None:
					bin_size = QueryEngine.NUM_BINS_3D
				query_string = QueryEngine.QUERY_TEMPLATES['3D_real_disc'].format(basis = basis, attr = bigquery_name, attr_name = gene_name, 
					x_min = self.df_ranges.at[basis + '_X', 'min_value'], x_denom = self.df_ranges.at[basis + '_X', 'max_value'] - self.df_ranges.at[basis + '_X', 'min_value'],
					y_min = self.df_ranges.at[basis + '_Y', 'min_value'], y_denom = self.df_ranges.at[basis + '_Y', 'max_value'] - self.df_ranges.at[basis + '_Y', 'min_value'],
					z_min = self.df_ranges.at[basis + '_Z', 'min_value'], z_denom = self.df_ranges.at[basis + '_Z', 'max_value'] - self.df_ranges.at[basis + '_Z', 'min_value'],
					bin_size = bin_size, dataset = self.dataset, table_id = table_id)
			num_rows = self.calc_num_rows(basis, bin_size)
		else:
			if basis != 'DIFFMAP':
				query_string = QueryEngine.QUERY_TEMPLATES['2D_real'].format(basis = basis, attr = bigquery_name, attr_name = gene_name, dataset = self.dataset, table_id = table_id)
			else:
				query_string = QueryEngine.QUERY_TEMPLATES['3D_real'].format(basis = basis, attr = bigquery_name, attr_name = gene_name, dataset = self.dataset, table_id = table_id)
			num_rows = self.nsample
		
		print(query_string)

		if self.is_parallel:
			results = self.query_parallel(query_string, num_rows)
		else:
			results = self.query_storage(query_string)

		return results

	def query_attribute(self, basis, attr, discretization = True, bin_size = None):
		results = None

		if attr not in self.attr2type:
			print("Query failed: {attr} does not exist!".format(attr = attr))
			return results

		qtype = 'cat' if self.attr2type[attr] == 'STRING' else 'real'

		# If sample size > MAX_SAMPLE_SIZE_NODISC, set discretization to True
		if not discretization and self.nsample > QueryEngine.MAX_SAMPLE_SIZE_NODISC:
			discretization = True

		query_string = ''
		num_rows = 0
		if discretization:
			if basis != 'DIFFMAP':
				if bin_size is None:
					bin_size = QueryEngine.NUM_BINS_2D
				query_string = QueryEngine.QUERY_TEMPLATES['2D_{qtype}_disc'.format(qtype = qtype)].format(basis = basis, attr = attr, attr_name = attr, 
					x_min = self.df_ranges.at[basis + '_X', 'min_value'], x_denom = self.df_ranges.at[basis + '_X', 'max_value'] - self.df_ranges.at[basis + '_X', 'min_value'],
					y_min = self.df_ranges.at[basis + '_Y', 'min_value'], y_denom = self.df_ranges.at[basis + '_Y', 'max_value'] - self.df_ranges.at[basis + '_Y', 'min_value'],
					bin_size = bin_size, dataset = self.dataset, table_id = 0)
			else:
				if bin_size is None:
					bin_size = QueryEngine.NUM_BINS_3D
				query_string = QueryEngine.QUERY_TEMPLATES['3D_{qtype}_disc'.format(qtype = qtype)].format(basis = basis, attr = attr, attr_name = attr, 
					x_min = self.df_ranges.at[basis + '_X', 'min_value'], x_denom = self.df_ranges.at[basis + '_X', 'max_value'] - self.df_ranges.at[basis + '_X', 'min_value'],
					y_min = self.df_ranges.at[basis + '_Y', 'min_value'], y_denom = self.df_ranges.at[basis + '_Y', 'max_value'] - self.df_ranges.at[basis + '_Y', 'min_value'],
					z_min = self.df_ranges.at[basis + '_Z', 'min_value'], z_denom = self.df_ranges.at[basis + '_Z', 'max_value'] - self.df_ranges.at[basis + '_Z', 'min_value'],
					bin_size = bin_size, dataset = self.dataset, table_id = 0)
			num_rows = self.calc_num_rows(basis, bin_size)
		else:
			if basis != 'DIFFMAP':
				query_string = QueryEngine.QUERY_TEMPLATES['2D_{qtype}'.format(qtype = qtype)].format(basis = basis, attr = attr, attr_name = attr, dataset = self.dataset, table_id = 0)
			else:
				query_string = QueryEngine.QUERY_TEMPLATES['3D_{qtype}'.format(qtype = qtype)].format(basis = basis, attr = attr, attr_name = attr, dataset = self.dataset, table_id = 0)
			num_rows = self.nsample
		
		print(query_string)

		if self.is_parallel:
			results = self.query_parallel(query_string, num_rows)
		else:
			results = self.query_storage(query_string)

		return results




