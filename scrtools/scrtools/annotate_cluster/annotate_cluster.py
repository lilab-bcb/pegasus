#!/usr/bin/env python

import numpy as np
import json
import pandas as pd
from sys import stdout


class CellType:
	def __init__(self, name):
		self.name = name
		self.score = self.avgp = 0.0
		self.weak_support = []
		self.strong_support = []
		self.subtypes = None

	def evaluate(self, obj, de_up, de_down, kwds, thre):
		self.score = self.avgp = 0.0
		self.weak_support = []
		self.strong_support = []

		nump = 0
		for marker_set in obj['markers']:
			numer = 0.0
			denom = len(marker_set['genes']) * 2.0
			if denom == 0.0:
				continue

			for marker in marker_set['genes']:
				sign = marker[-1]
				gsym = marker[:-1]

				if sign == '+':
					if gsym in de_up.index:
						fc = de_up.at[gsym, kwds['fc']]
						expr = de_up.at[gsym, kwds['expr']]
						self.avgp += expr
						nump += 1

						if fc >= thre:
							numer += 2.0
							self.strong_support.append((marker, "{0:.2f}".format(expr)))
						else:
							numer += 1.0 + (fc - 1.0) / (thre - 1.0)
							self.weak_support.append((marker, "{0:.2f}".format(expr)))
				else:
					assert sign == '-'
					if gsym not in de_up.index:
						if gsym in de_down.index:
							fc = (1.0 / de_down.at[gsym, kwds['fc']]) if de_down.at[gsym, kwds['fc']] > 0.0 else np.inf
							expr = de_down.at[gsym, kwds['expr']]
							if fc >= thre:
								numer += 2.0
								self.strong_support.append((marker, "{0:.2f}".format(expr)))
							else:
								numer += 1.0 + (fc - 1.0) / (thre - 1.0)
								self.weak_support.append((marker, "{0:.2f}".format(expr)))
						else:
							numer += 1.0
							self.weak_support.append((marker, "N/A"))
			
			self.score += numer / denom * marker_set['weight']

		self.score = self.score / obj['denominator'] if obj['denominator'] > 0.0 else 0.0
		if nump > 0:
			self.avgp /= nump

	def tostring(self):
		res = "name: {0}; score: {1:.2f}; avgp: {2:.2f}".format(self.name, self.score, self.avgp)
		if len(self.strong_support) > 0:
			res += "; strong support: {0}".format(",".join(["({0},{1})".format(x[0], x[1]) for x in self.strong_support]))
		if len(self.weak_support) > 0:
			res += "; weak support: {0}".format(",".join(["({0},{1})".format(x[0], x[1]) for x in self.weak_support]))
		
		return res


class Annotator:
	def __init__(self, json_file, genes):
		with open(json_file) as fin:
			self.object = json.load(fin)
		self.recalibrate(self.object, genes)	

	def recalibrate(self, obj, genes):
		for celltype in obj['cell_types']:
			denom = 0.0
			for marker_set in celltype['markers']:
				markers = marker_set['genes']
				s = len(markers)
				marker_set['genes'] = [x for x in markers if x[:-1] in genes]
				new_s = len(marker_set['genes'])
				marker_set['weight'] = marker_set['weight'] / s * new_s
				denom += marker_set['weight']
			celltype['denominator'] = denom
			sub_obj = celltype.get('subtypes', None)
			if sub_obj is not None:
				self.recalibrate(sub_obj, genes)
	
	def evaluate(self, de_up, de_down, kwds, thre = 1.5, obj = None):
		if obj is None:
			obj = self.object

		results = []
		for celltype in obj['cell_types']:
			ct = CellType(celltype['name'])
			ct.evaluate(celltype, de_up, de_down, kwds, thre)
			if ct.score > 0.5:
				sub_obj = celltype.get('subtypes', None)
				if sub_obj is not None:
					ct.subtypes = self.evaluate(de_up, de_down, kwds, thre, sub_obj)
			results.append(ct)

		results.sort(key = lambda x: x.score, reverse = True)

		return results
	
	def report(self, fout, results, thre, space = 4):
		for ct in results:
			if ct.score > thre: 
				fout.write(' ' * space + ct.tostring() + '\n')
				if ct.subtypes is not None:
					self.report(fout, ct.subtypes, 0.5, space + 4)


def annotate_clusters(adata, test, json_file, thre = 0.5, fout = stdout):
	anno = Annotator(json_file, adata.var_names)
	if test == "fisher":
		kwds = {'fc' : 'fold_change', 'expr' : 'percentage'}
	else:
		kwds = {'fc' : 'log_fold_change', 'expr' : 'mean_log_expression'}
	size = adata.obs['louvain_groups'].cat.categories.size
	for i in range(size):
		clust_str = "de_{test}_{clust}".format(test = test, clust = i)
		de_up = pd.DataFrame(data = adata.uns[clust_str + "_up_stats"], index = pd.Index(data = adata.uns[clust_str + "_up_genes"], name = "gene"))
		de_down = pd.DataFrame(data = adata.uns[clust_str + "_down_stats"], index = pd.Index(data = adata.uns[clust_str + "_down_genes"], name = "gene"))
		results = anno.evaluate(de_up, de_down, kwds)
		fout.write("Cluster {0}:\n".format(i + 1))
		anno.report(fout, results, thre)
