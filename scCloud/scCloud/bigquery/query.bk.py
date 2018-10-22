import numpy as np
import pandas as pd
import uuid
import time
import threading
from google.cloud import bigquery

class QueryEngine:
	def __init__(self, project_id):
		# if unique_id is None:
		# 	unique_id = uuid.uuid4().hex

		self.project_id = project_id
		self.client = bigquery.Client(project = project_id)

		# self.storage_client = storage.Client(project = project_id)
		# try:
		# 	self.bucket = self.storage_client.get_bucket(self.unique_id)
		# except google.cloud.exceptions.NotFound:
		# 	self.bucket = self.storage_client.create_bucket(self.unique_id)

		# self.dataset_ref = self.client.dataset(self.unique_id)
		# self.table_ref = self.dataset_ref.table(self.table_name)

		# dataset = bigquery.Dataset(self.dataset_ref)
		# dataset.location = 'US'
		# try:
		# 	dataset = self.client.create_dataset(dataset)
		# except Exception:
		# 	None

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		print(exc_type)
		print(exc_value)
		print(traceback)
		# try:
		# 	self.client.delete_dataset(self.dataset_ref, delete_contents = True)
		# except Exception:
		# 	None

	# def query(query):
	# 	# set up job configuration
	# 	query_config = bigquery.QueryJobConfig()
	# 	query_config.allow_large_results = True
	# 	query_config.destination = self.table_ref
	# 	query_config.write_disposition = 'WRITE_TRUNCATE'

	# 	# sent query to BigQuery
	# 	query_job = self.client.query(query = query, location = 'US', job_config = query_config)
	# 	query_job.result()

	# 	# extract to google bucket
	# 	extract_job = self.client.extract_table(self.table_ref, 'gs://{bucket}/vis_results.csv'.format(bucket = self.bucket.id))
	# 	extract_job.result()

	def collect_data(self, thread_no, num_base, num_left, df_arr, table_ref):

		start_index = num_base * thread_no + min(thread_no, num_left)
		max_results = num_base + (thread_no < num_left)

		client = bigquery.Client(project = self.project_id)
		result_table = client.get_table(client.dataset(table_ref.dataset_id).table(table_ref.table_id))
		rows = client.list_rows(result_table, start_index = start_index, max_results = max_results)

		print("~ Begint thread {}".format(thread_no))
		
		start_time = time.time()
		for i, page in enumerate(rows.pages):
			print("Thread {}: {}, {}, {}".format(thread_no, i, rows.page_number, page.num_items))
		# df_arr[thread_no] = rows.to_dataframe()

		end_time = time.time()
		print("Thread = {}, Time = {:.2f}.".format(thread_no, end_time - start_time))

	def query(self, query, num_threads = 10, page_size = 100000):
		start = time.time()
		# sent query to BigQuery
		query_job = self.client.query(query = query, location = 'US')
		query_job.result()
		table_ref = query_job.destination
		result_table = self.client.get_table(table_ref)

		end = time.time()
		print("Time spent for querying = {:.2f}s.".format(end - start))
		start = end

		num_base = result_table.num_rows // num_threads
		num_left = result_table.num_rows % num_threads

		threads = [None] * num_threads
		df_arr = [None] * num_threads
		for i in range(num_threads):
			t = threading.Thread(target=self.collect_data, args=(i, num_base, num_left, df_arr, table_ref))
			threads[i] = t
			t.start()

		for i in range(num_threads):
			threads[i].join()

		end = time.time()
		print("Time spent for loading = {:.2f}s".format(end - start))
		start = time.time()

		df = None
		# df = pd.concat(df_arr)

		end = time.time()
		print("Time spent = {:.2f}s.".format(end - start))

		return df
