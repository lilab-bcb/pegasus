import time
import numpy as np
import pandas as pd
from collections import defaultdict
import xlsxwriter

from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
# from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

from . import read_input

def find_markers(data, label_attr, n_jobs = 1, min_gain = 1.0, random_state = 0, remove_ribo = False):
	if remove_ribo:
		data = data[:,np.vectorize(lambda x: not x.startswith('RPL') and not x.startswith('RPS'))(data.var_names)]

	X_train, X_test, y_train, y_test = train_test_split(data.X, data.obs[label_attr], test_size = 0.1, random_state = random_state, stratify = data.obs[label_attr])

	# start = time.time()
	# xgb = XGBClassifier(n_jobs = n_jobs, n_gpus = 0)
	# xgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], eval_metric = 'merror')
	# # print(xgb.evals_result())
	# end = time.time()
	# print("XGBoost used {:.2f}s to train.".format(end - start))

	start = time.time()
	lgb = LGBMClassifier(n_jobs = n_jobs, metric = 'multi_error', importance_type = 'gain')
	lgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], early_stopping_rounds = 1)
	end = time.time()
	print("LightGBM used {:.2f}s to train.".format(end - start))

	ntot = (lgb.feature_importances_ >= min_gain).sum()
	ords = np.argsort(lgb.feature_importances_)[::-1][:ntot]

	ncat = data.obs[label_attr].cat.categories.size
	log_exprs = ['mean_log_expression_{}'.format(i + 1) for i in range(ncat)]
	titles = [('down', 'down_gain'), ('weak', 'weak_gain'), ('strong', 'strong_gain')]
	markers = defaultdict(lambda: defaultdict(list))

	kmeans = KMeans(n_clusters = 3, random_state = random_state)
	for gene_id in ords:
		gene_symbol = data.var_names[gene_id]
		mydat = data.var.loc[gene_symbol, log_exprs].values.reshape(-1, 1)
		kmeans.fit(mydat)
		kmeans_label_mode = pd.Series(kmeans.labels_).mode()[0]
		for i, kmeans_label in enumerate(np.argsort(kmeans.cluster_centers_[:,0])):
			if kmeans_label != kmeans_label_mode:
				for clust_label in (kmeans.labels_ == kmeans_label).nonzero()[0]:
					markers[clust_label][titles[i][0]].append(gene_symbol)
					markers[clust_label][titles[i][1]].append('{:.2f}'.format(lgb.feature_importances_[gene_id]))

	end = time.time()
	print("find_markers took {:.2f}s to finish.".format(end - start))

	return markers



def run_find_markers(input_h5ad_file, output_file, label_attr = 'louvain_labels', n_jobs = 1, min_gain = 1.0, random_state = 0, remove_ribo = False):
	data = read_input(input_h5ad_file, mode = 'a')
	markers = find_markers(data, label_attr, n_jobs = n_jobs, min_gain = min_gain, random_state = random_state, remove_ribo = remove_ribo)
	
	nclust = len(markers)
	keywords = [('strong', 'strong_gain'), ('weak', 'weak_gain'), ('down', 'down_gain')]

	writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
	
	for i in range(nclust):
		sizes = []
		for keyword in keywords:
			sizes.append(len(markers[i][keyword[0]]))
		
		arr = np.zeros((max(sizes), 8), dtype = object)
		arr[:] = ''
		
		for j in range(3):
			arr[0:sizes[j], j * 3] = markers[i][keywords[j][0]]
			arr[0:sizes[j], j * 3 + 1] = markers[i][keywords[j][1]]
		
		df = pd.DataFrame(data = arr, columns = ['strongly up-regulated', 'gain', '', 'weakly up-regulated', 'gain', '', 'down-regulated', 'gain'])
		df.to_excel(writer, sheet_name = "{}".format(i + 1), index = False)

	writer.save()	
