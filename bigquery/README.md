# Bigquery codes in Node.js

## Usage:

```
const QueryEngine = require('query.js');
const qe = new QueryEngine(projectId, datasetId, [discretization = true]);
res_obj = qe.query('TSNE', 'CD8A', true); // basis, attr, is_gene
res_obj.promise.then(() => {
    do something with res_obj.results.
  });
```

## res_obj.results

This is an object. Bases and attr/gene are its properties. Each property is an array of values that people can fit into plotly.

