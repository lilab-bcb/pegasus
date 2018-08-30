const BigQuery = require('@google-cloud/bigquery');
const Storage = require('@google-cloud/storage');
const fillTemplate = require('es6-dynamic-template');
const fs = require('fs');

const bucketName = 'a6fc7179-6fe5-420b-b895-89a73d7250ee'; // default bucket name

class QueryEngine {
	constructor(project_id, dataset, discretization = true, location = 'US') {
		this.project_id = project_id;
		this.dataset = dataset;
		this.discretization = discretization;
		this.location = location;

		this.bigquery = new BigQuery({projectId: this.project_id});
		this.storage = new Storage({projectId: this.project_id});
		
		var qe_obj = this;
		var start_time = Date.now();
		var query, options, promise;

		this.promises = []

		promise = this.storage.bucket(bucketName)
			.exists()
			.then(data => {
				if (!data[0]) {
					qe_obj.storage
						  .createBucket(bucketName)
						  .then(() => {
							console.log(`Bucket ${bucketName} created.`);	
						  })
						  .catch(err => {
						  	console.error('ERROR:', err);
						  });
				}
			});
		this.promises.push(promise);

		query = fillTemplate(QueryEngine.QUERY_TEMPLATES.ranges, {dataset : this.dataset});
		options = {query: query, useLegacySql: false};
		promise = this.bigquery.query(options)
					 .then(results => {
							qe_obj.ranges = QueryEngine.rows_to_obj(results[0], 'name', ['min_value', 'max_value']);
							console.log('Time spent for querying ranges = ' + (Date.now() - start_time) / 1000.0 + ' s.')
  					  })
					 .catch(err => {
    						console.error('ERROR:', err);
  					  });
		this.promises.push(promise);


		query = fillTemplate(QueryEngine.QUERY_TEMPLATES.gene_names, {dataset : this.dataset});
		options = {query: query, useLegacySql: false};
		promise = this.bigquery.query(options)
					 .then(results => {
							qe_obj.gene_names = QueryEngine.rows_to_obj(results[0], 'gene', ['bigquery_name', 'table_number']);
							console.log('Time spent for querying gene_names = ' + (Date.now() - start_time) / 1000.0 + ' s.')
  					  })
					 .catch(err => {
    						console.error('ERROR:', err);
  					  });
		this.promises.push(promise);

		promise = this.bigquery
						.dataset(this.dataset)
						.table('expression_0')
						.getMetadata()
						.then(data => {
							qe_obj.nsample = data[0].numRows;
						})
						.catch(err => {
							console.error('ERROR:', err);
						});
		this.promises.push(promise);		
	}

	generate_query(basis, attr, is_gene) {
		var query;

		if (is_gene) {
			if (basis !== 'DIFFMAP') {
				if (this.discretization) {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D2_real_disc,
										{
											basis : basis, 
											attr : this.gene_names[attr].bigquery_name,
								 			attr_name : attr,
								 			x_min : this.ranges[basis + '_X'].min_value,
								 			x_denom : this.ranges[basis + '_X'].max_value - this.ranges[basis + '_Y'].min_value,
								 			y_min : this.ranges[basis + '_Y'].min_value,
								 			y_denom : this.ranges[basis + '_Y'].max_value - this.ranges[basis + '_Y'].min_value,
								 			bin_size : QueryEngine.NUM_BINS_D2,
								 			dataset : this.dataset,
								 			table_id : this.gene_names[attr].table_number
										});
				} else {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D2_real,
										{
											basis : basis, 
											attr : this.gene_names[attr].bigquery_name,
								 			attr_name : attr,
								 			dataset : this.dataset,
								 			table_id : this.gene_names[attr].table_number
										});
				}
			} else {
				if (this.discretization) {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D3_real_disc,
										{
											basis : basis, 
											attr : this.gene_names[attr].bigquery_name,
								 			attr_name : attr,
								 			x_min : this.ranges[basis + '_X'].min_value,
								 			x_denom : this.ranges[basis + '_X'].max_value - this.ranges[basis + '_Y'].min_value,
								 			y_min : this.ranges[basis + '_Y'].min_value,
								 			y_denom : this.ranges[basis + '_Y'].max_value - this.ranges[basis + '_Y'].min_value,
								 			z_min : this.ranges[basis + '_Z'].min_value,
								 			z_denom : this.ranges[basis + '_Z'].max_value - this.ranges[basis + '_Z'].min_value,
								 			bin_size : QueryEngine.NUM_BINS_D3,
								 			dataset : this.dataset,
								 			table_id : this.gene_names[attr].table_number
										});
				} else {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D3_real,
										{
											basis : basis, 
											attr : this.gene_names[attr].bigquery_name,
								 			attr_name : attr,
								 			dataset : this.dataset,
								 			table_id : this.gene_names[attr].table_number
										});
				}
			}
		} else {
			if (basis !== 'DIFFMAP') {
				if (this.discretization) {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D2_cat_disc,
										{
											basis : basis, 
											attr : attr,
								 			attr_name : attr,
								 			x_min : this.ranges[basis + '_X'].min_value,
								 			x_denom : this.ranges[basis + '_X'].max_value - this.ranges[basis + '_Y'].min_value,
								 			y_min : this.ranges[basis + '_Y'].min_value,
								 			y_denom : this.ranges[basis + '_Y'].max_value - this.ranges[basis + '_Y'].min_value,
								 			bin_size : QueryEngine.NUM_BINS_D2,
								 			dataset : this.dataset,
								 			table_id : 0
										});
				} else {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D2_cat,
										{
											basis : basis, 
											attr : attr,
								 			attr_name : attr,
								 			dataset : this.dataset,
								 			table_id : 0
										});
				}
			} else {
				if (this.discretization) {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D3_cat_disc,
										{
											basis : basis, 
											attr : attr,
								 			attr_name : attr,
								 			x_min : this.ranges[basis + '_X'].min_value,
								 			x_denom : this.ranges[basis + '_X'].max_value - this.ranges[basis + '_Y'].min_value,
								 			y_min : this.ranges[basis + '_Y'].min_value,
								 			y_denom : this.ranges[basis + '_Y'].max_value - this.ranges[basis + '_Y'].min_value,
								 			z_min : this.ranges[basis + '_Z'].min_value,
								 			z_denom : this.ranges[basis + '_Z'].max_value - this.ranges[basis + '_Z'].min_value,
								 			bin_size : QueryEngine.NUM_BINS_D3,
								 			dataset : this.dataset,
								 			table_id : 0
										});
				} else {
					query = fillTemplate(QueryEngine.QUERY_TEMPLATES.D3_cat,
										{
											basis : basis, 
											attr : attr,
								 			attr_name : attr,
								 			dataset : this.dataset,
								 			table_id : 0
										});
				}
			}
		}

		return query;
	}

	query(basis, attr, is_gene) {
		var qe_obj = this;
		var res_obj = {};
		var dest; // temp table storing query results

		res_obj.promise = Promise.all(this.promises).then(values => {
			const query = qe_obj.generate_query(basis, attr, is_gene);
			const options = {
				query: query,
				useLegacySql: false
			};

			var start_time = Date.now();
			var pms;

			if (qe_obj.nsample > QueryEngine.MAX_SAMPLE_SIZE_DIRECT) {
				let csv_file = qe_obj.storage.bucket(bucketName).file('temp_vis.csv');

				pms = qe_obj.bigquery
		  					.createQueryJob(options)
							.then(results => {
	  							var job = results[0];
	  							dest = results[1].configuration.query.destinationTable;
	  							return job.promise();
		  					})
		  					.then(() => {
		  						return qe_obj.bigquery
		  										.dataset(dest.datasetId)
		  										.table(dest.tableId)
		  										.extract(csv_file)
				  								.then(results => {
				  									var job = results[0];
				  									return job;
				  								})
				  								.then(() => {
				  									return csv_file.download().then(data => {
				  										res_obj.results = QueryEngine.parse_csv_by_column(data[0], basis, attr, is_gene);
				  										console.log(`Time spent for querying ${attr} = ` + (Date.now() - start_time) / 1000.0 + ' s.');
				  									})
				  								});
				  									
		  					});
		  	} else {
		  		pms = qe_obj.bigquery
			  					.query(options)
								.then(results => {
		  							res_obj.results = QueryEngine.rows_by_column(results[0], basis, attr);
		  							console.log(`Time spent for querying ${attr} = ` + (Date.now() - start_time2) / 1000.0 + ' s.');
			  					});
		  	}

		  	return pms;
		})
		.catch(err => {
			console.error('ERROR:', err);
		});

		return res_obj;
	}


	static rows_to_obj(rows, prop_key, value_keys) {
		var obj = {}
		for (var i = 0; i < rows.length; ++i) {
			let prop = rows[i][prop_key]
			obj[prop] = {}
			value_keys.forEach(function(item) {obj[prop][item] = rows[i][item]})
		}
		return obj
	}

	static rows_by_column(rows, basis, attr) {
		var results = {};
		var columns = [];

		const arow = rows[0];
		for (let key in arow)
			if (arow.hasOwnProperty(key) && (key === attr || key.indexOf(basis) === 0)) {
				columns.push(key);
				results[key] = [];
			}
		rows.forEach(function(row) {
			columns.forEach(function(key) {
				results[key].push(row[key]);
			});
		});

		return results;
	}

	static parse_csv_by_column(data, basis, attr, is_gene) {
		var results = {};
		const lines = data.toString().split(/\r\n|\n/);

		var posvec = [];
		var columns = [];
		var converts = [];
		lines[0].split(',').forEach((item, i) => {
			if (item.indexOf(basis) === 0) {
				results[item] = [];
				posvec.push(i);
				columns.push(item);
				converts.push(true);
			} else if (item === attr) {
				results[item] = [];
				posvec.push(i);
				columns.push(item);
				converts.push(is_gene ? true : false);				
			}
		});

		lines.slice(1).forEach(line => {
			let fields = line.split(',');
			posvec.forEach((pos, i) => {
				results[columns[i]].push(converts[i] ? parseFloat(fields[pos]) : fields[pos]);
			})
		})

		return results;
	}
}

QueryEngine.NUM_BINS_D2 = 1024; // default D2 discretization granule
QueryEngine.NUM_BINS_D3 = 512; // default D3 discretization granule
QueryEngine.MAX_SAMPLE_SIZE_DIRECT = 100000; // maximum sample size to use direct query

QueryEngine.QUERY_TEMPLATES = {
	gene_names : `SELECT * FROM \\\`\${dataset}.gene_names\\\``,
	ranges : `SELECT * FROM \\\`\${dataset}.ranges\\\``,

	D2_real : `SELECT ROUND(\${basis}_X, 5) AS \${basis}_X, ROUND(\${basis}_Y, 5) AS \${basis}_Y, ROUND(\${attr}, 2) AS \${attr_name} FROM \\\`\${dataset}.expression_\${table_id}\\\``,
	D2_cat : `SELECT ROUND(\${basis}_X, 5) AS \${basis}_X, ROUND(\${basis}_Y, 5) AS \${basis}_Y, \${attr} AS \${attr_name} FROM \\\`\${dataset}.expression_\${table_id}\\\``,
	D3_real : `SELECT ROUND(\${basis}_X, 5) AS \${basis}_X, ROUND(\${basis}_Y, 5) AS \${basis}_Y, ROUND(\${basis}_Z, 5) AS \${basis}_Z, ROUND(\${attr}, 2) AS \${attr_name} FROM \\\`\${dataset}.expression_\${table_id}\\\``,
	D3_cat : `SELECT ROUND(\${basis}_X, 5) AS \${basis}_X, ROUND(\${basis}_Y, 5) AS \${basis}_Y, ROUND(\${basis}_Z, 5) AS \${basis}_Z, \${attr} AS \${attr_name} FROM \\\`\${dataset}.expression_\${table_id}\\\``,
	
	D2_real_disc : `SELECT ROUND(AVG(\${basis}_X), 5) AS \${basis}_X,
							 ROUND(AVG(\${basis}_Y), 5) AS \${basis}_Y,
							 ROUND(AVG(\${attr}), 2) AS \${attr_name},
							 CAST(FLOOR((\${basis}_X - \${x_min}) / \${x_denom} * \${bin_size}) AS INT64) * \${bin_size} + CAST(FLOOR((\${basis}_Y - \${y_min}) / \${y_denom} * \${bin_size}) AS INT64) AS MYID
					  FROM \\\`\${dataset}.expression_\${table_id}\\\`
					  GROUP BY MYID`,
	D2_cat_disc : `SELECT ROUND(AVG(\${basis}_X), 5) AS \${basis}_X,
					 		ROUND(AVG(\${basis}_Y), 5) AS \${basis}_Y,
							APPROX_TOP_COUNT(\${attr}, 1)[OFFSET(0)].value AS \${attr_name},
							CAST(FLOOR((\${basis}_X - \${x_min}) / \${x_denom} * \${bin_size}) AS INT64) * \${bin_size} + CAST(FLOOR((\${basis}_Y - \${y_min}) / \${y_denom} * \${bin_size}) AS INT64) AS MYID
					 FROM \\\`\${dataset}.expression_\${table_id}\\\` 
					 GROUP BY MYID`,

	D3_real_disc : `SELECT ROUND(AVG(\${basis}_X), 5) AS \${basis}_X,
							 ROUND(AVG(\${basis}_Y), 5) AS \${basis}_Y,
							 ROUND(AVG(\${basis}_Z), 5) AS \${basis}_Z,
							 ROUND(AVG(\${attr}), 2) AS \${attr_name},
							 CAST(FLOOR((\${basis}_X - \${x_min}) / \${x_denom} * \${bin_size}) AS INT64) * \${bin_size} * \${bin_size} + CAST(FLOOR((\${basis}_Y - \${y_min}) / \${y_denom} * \${bin_size}) AS INT64) * \${bin_size} + CAST(FLOOR((\${basis}_Z - \${z_min}) / \${z_denom} * \${bin_size}) AS INT64) AS MYID
					  FROM \\\`\${dataset}.expression_\${table_id}\\\`
					  GROUP BY MYID`,
	D3_cat_disc : `SELECT ROUND(AVG(\${basis}_X), 5) AS \${basis}_X,
							ROUND(AVG(\${basis}_Y), 5) AS \${basis}_Y,
							ROUND(AVG(\${basis}_Z), 5) AS \${basis}_Z,
							APPROX_TOP_COUNT(\${attr}, 1)[OFFSET(0)].value AS \${attr_name},
							CAST(FLOOR((\${basis}_X - \${x_min}) / \${x_denom} * \${bin_size}) AS INT64) * \${bin_size} * \${bin_size} + CAST(FLOOR((\${basis}_Y - \${y_min}) / \${y_denom} * \${bin_size}) AS INT64) * \${bin_size} + CAST(FLOOR((\${basis}_Z - \${z_min}) / \${z_denom} * \${bin_size}) AS INT64) AS MYID
					 FROM \\\`\${dataset}.expression_\${table_id}\\\`
					 GROUP BY MYID`
}

module.exports = QueryEngine;
