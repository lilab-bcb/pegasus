import time
import numpy as np
import pandas as pd
import xlsxwriter

from xgboost import XGBClassifier


def find_markers_with_xgboost(data, label_attr, output_name, n_jobs = 1, top_n = 500, alpha = 0.05)
	start = time.time()
	xgb = XGBClassifier(n_jobs = n_jobs)
	xgb.fit(data.X, data.obs[label_attr])
	end = time.time()
	print("Spent {:.2f}s in training XGBoost model.".format(end - start))

	xgb.save_model(output_name + '.xgb.model')
	y_pred = xgb.predict(data.X)
	print("Accuracy = {:.2%}.".format((y_pred == data.obs[label_attr]).sum() / data.shape[0]))

	top_n = min(top_n, (xgb.feature_importances_ > 0).sum())
	genes = data.var_names.values[np.argsort(xgb.feature_importances_)[::-1][:top_n]]

	writer = pd.ExcelWriter(output_name + 'xgb.markers.xlsx' engine='xlsxwriter')

	var_df = data.var.loc[genes]
	cols = ["percentage", "percentage_other", "percentage_fold_change", "mean_log_expression", "log_fold_change", "WAD_score", "auc", "predpower"]	
	for clust_id in data.obs[label_attr]:
		idx = var_df["t_pval_{0}".format(clust_id)] <= alpha
		idx_up = idx & (var_df["WAD_score_{0}".format(clust_id)] > 0.0)
		idx_down = idx & (var_df["WAD_score_{0}".format(clust_id)] < 0.0)
		assert idx_up.sum() + idx_down.sum() == idx.sum()

		col_names = ["{0}_{1}".format(x, clust_id) for x in cols]

		df_up = pd.DataFrame(var_df.loc[idx_up.values, col_names])
		df_up.rename(columns = lambda x: '_'.join(x.split('_')[:-1]), inplace = True)
		df_up.to_excel(writer, sheet_name = "{0} up".format(cluster_id))
								
		df_down = pd.DataFrame(var_df.loc[idx_down.values, col_names])
		df_down.rename(columns = lambda x: '_'.join(x.split('_')[:-1]), inplace = True)
		df_down.to_excel(writer, sheet_name = "{0} down".format(cluster_id))

	writer.save()
