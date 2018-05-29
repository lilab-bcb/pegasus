import numpy as np
import json
import pandas as pd
from sys import stdout
from natsort import natsorted

class CellType:
	def __init__(self, name, ignoreNA = False):
		self.name = name
		self.score = self.avgp = 0.0
		self.weak_support = []
		self.strong_support = []
		self.subtypes = None
		self.ignoreNA = ignoreNA

	def evaluate(self, obj, de_up, de_down, thre):
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
						fc = de_up.at[gsym, 'fc']
						expr = de_up.at[gsym, 'expr']
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
							fc = (1.0 / de_down.at[gsym, 'fc']) if de_down.at[gsym, 'fc'] > 0.0 else np.inf
							expr = de_down.at[gsym, 'expr']
							if fc >= thre:
								numer += 2.0
								self.strong_support.append((marker, "{0:.2f}".format(expr)))
							else:
								numer += 1.0 + (fc - 1.0) / (thre - 1.0)
								self.weak_support.append((marker, "{0:.2f}".format(expr)))
						elif not self.ignoreNA:
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
	
	def evaluate(self, de_up, de_down, thre = 1.5, obj = None, ignoreNA = False):
		if obj is None:
			obj = self.object

		results = []
		for celltype in obj['cell_types']:
			ct = CellType(celltype['name'], ignoreNA)
			ct.evaluate(celltype, de_up, de_down, thre)
			if ct.score > 0.5:
				sub_obj = celltype.get('subtypes', None)
				if sub_obj is not None:
					ct.subtypes = self.evaluate(de_up, de_down, thre, sub_obj, ignoreNA = ignoreNA)
			results.append(ct)

		results.sort(key = lambda x: x.score, reverse = True)

		return results
	
	def report(self, fout, results, thre, space = 4):
		for ct in results:
			if ct.score > thre: 
				fout.write(' ' * space + ct.tostring() + '\n')
				if ct.subtypes is not None:
					self.report(fout, ct.subtypes, 0.5, space + 4)



def annotate_clusters(data, json_file, thre, fout = stdout, ignoreNA = False):
	anno = Annotator(json_file, data.var_names)
	
	clusts = natsorted([x[10:] for x in data.var.columns if x.startswith("WAD_score_")])
	tests = [x for x in ['t', 'fisher', 'mwu'] if "{0}_qval_{1}".format(x, clusts[0]) in data.var.columns]
	for clust_id in clusts:
		idx = data.var["{0}_qval_{1}".format(tests[0], clust_id)] <= 0.05
		for test in tests[1:]:
			idx = idx & (data.var["{0}_qval_{1}".format(test, clust_id)] <= 0.05)

		idx_up = idx & (data.var["WAD_score_{0}".format(clust_id)] > 0.0)
		idx_down = idx & (data.var["WAD_score_{0}".format(clust_id)] < 0.0)
		assert idx_up.sum() + idx_down.sum() == idx.sum()

		cols = ["{0}_{1}".format(x, clust_id) for x in ["percentage_fold_change", "percentage"]]
		de_up = pd.DataFrame(data.var.loc[idx_up.values, cols])
		de_up.rename(columns = {cols[0]: "fc", cols[1]: "expr"}, inplace = True)
		de_down = pd.DataFrame(data.var.loc[idx_down.values, cols])
		de_down.rename(columns = {cols[0]: "fc", cols[1]: "expr"}, inplace = True)
		results = anno.evaluate(de_up, de_down, ignoreNA = ignoreNA)
		fout.write("Cluster {0}:\n".format(clust_id))
		anno.report(fout, results, thre)
