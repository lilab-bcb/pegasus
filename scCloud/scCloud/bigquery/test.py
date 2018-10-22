#!/usr/bin/env python

import pandas_gbq as gbq




idx_df = data.obs_names.isin(adt.obs_names)
if idx_df.sum() < data.shape[0]:
	nzero = data.shape[0] - idx_df.sum()
	print("Warning: {} cells do not have ADTs, percentage = {:.2f}%.".format(nzero, nzero * 100.0 / data.shape[0]))

adt_small = adt[data.obs_names[idx_df],].X.toarray()
value = np.median(adt_small.sum(axis = 1))
tmp = adt[adt.obs_names.difference(data.obs_names),]
idx = tmp.X.sum(axis = 1).A1 < value
pvec = tmp.X[idx, ].sum(axis = 0).A1
pvec /= pvec.sum()

data.obsm['raw_probs'] = np.zeros((data.shape[0], adt.shape[1] + 1))
data.obsm['raw_probs'][:, adt.shape[1]] = 1.0




data3 = scrtools.tools.read_input('mouse_cortex_hashtag2_raw_h5.h5', 'mm10-1.2.0_premrna')
scrtools.tools.update_var_names(data3, 'mm10-1.2.0_premrna')
scrtools.tools.filter_data(data3, mito_prefix = 'mt-', min_genes = 200)
scrtools.tools.log_norm(data3, 1e5)
adt3 = scrtools.tools.read_input('mouse_cortex_hashtag2_ADTs.csv')
scc.demultiplex(data3, adt3)


start = time.time()
tmp2 = tmp[:, 0:9000].toarray()
end = time.time()
print(end - start)

from google.cloud import bigquery
client = bigquery.Client('studious-sign-138705')
dataset_ref = client.dataset('mydata')
query = ('SELECT * FROM `mydata.expr0` LIMIT 10')
query_job = client.query(query, location='US')
