import numpy as np
import pandas as pd
import json
import logging

from sys import stdout
from natsort import natsorted
from typing import List, Dict, Union

logger = logging.getLogger('sccloud')

class CellType:
    def __init__(self, name: str, ignore_nonde: bool = False):
        self.name = name
        self.score = self.avgp = 0.0
        self.weak_support = []
        self.strong_support = []
        self.subtypes = None
        self.ignore_nonde = ignore_nonde

    def evaluate(self, obj: 'json object', de_up: pd.DataFrame, de_down: pd.DataFrame, thre: float):
        """ Calculate score for matching a cluster with a putative cell type.
        """
        self.score = self.avgp = 0.0
        self.weak_support = []
        self.strong_support = []

        nump = 0
        for marker_set in obj["markers"]:
            numer = 0.0
            denom = len(marker_set["genes"]) * 2.0
            if denom == 0.0:
                continue

            for marker in marker_set["genes"]:
                sign = marker[-1]
                gsym = marker[:-1]

                if sign == "+":
                    if gsym in de_up.index:
                        fc = de_up.at[gsym, "fc"]
                        percent = de_up.at[gsym, "percent"]
                        self.avgp += percent
                        nump += 1

                        if fc >= thre:
                            numer += 2.0
                            self.strong_support.append((marker, "{0:.2f}%".format(percent)))
                        else:
                            numer += 1.0 + (fc - 1.0) / (thre - 1.0)
                            self.weak_support.append((marker, "{0:.2f}%".format(percent)))
                else:
                    assert sign == "-"
                    if gsym not in de_up.index:
                        if gsym in de_down.index:
                            fc = (
                                (1.0 / de_down.at[gsym, "fc"])
                                if de_down.at[gsym, "fc"] > 0.0
                                else np.inf
                            )
                            percent = de_down.at[gsym, "percent"]
                            if fc >= thre:
                                numer += 2.0
                                self.strong_support.append(
                                    (marker, "{0:.2f}%".format(percent))
                                )
                            else:
                                numer += 1.0 + (fc - 1.0) / (thre - 1.0)
                                self.weak_support.append(
                                    (marker, "{0:.2f}%".format(percent))
                                )
                        elif not self.ignore_nonde:
                            numer += 1.0
                            self.weak_support.append((marker, "N/A"))

            self.score += numer / denom * marker_set["weight"]

        self.score = (
            self.score / obj["denominator"] if obj["denominator"] > 0.0 else 0.0
        )
        if nump > 0:
            self.avgp /= nump

    def tostring(self):
        res = "name: {0}; score: {1:.2f}; average marker percentage: {2:.2f}%".format(
            self.name, self.score, self.avgp
        )
        if len(self.strong_support) > 0:
            res += "; strong support: {0}".format(
                ",".join(["({0},{1})".format(x[0], x[1]) for x in self.strong_support])
            )
        if len(self.weak_support) > 0:
            res += "; weak support: {0}".format(
                ",".join(["({0},{1})".format(x[0], x[1]) for x in self.weak_support])
            )

        return res


class Annotator:
    def __init__(self, marker_file: Union[str, Dict], genes: List[str]) -> None:
        if type(marker_file) != dict:
            with open(marker_file) as fin:
                self.object = json.load(fin)
        else:
            self.object = marker_file
        self.recalibrate(self.object, genes)

    def recalibrate(self, obj: 'json object', genes: List[str]) -> None:
        """ Remove markers that are not expressed (not in genes) and calculate partial weights for existing genes.
        """
        for celltype in obj["cell_types"]:
            denom = 0.0
            for marker_set in celltype["markers"]:
                markers = marker_set["genes"]
                s = len(markers)
                marker_set["genes"] = [x for x in markers if x[:-1] in genes]
                new_s = len(marker_set["genes"])
                marker_set["weight"] = marker_set["weight"] / s * new_s
                denom += marker_set["weight"]
            celltype["denominator"] = denom
            sub_obj = celltype.get("subtypes", None)
            if sub_obj is not None:
                self.recalibrate(sub_obj, genes)

    def evaluate(self, de_up: pd.DataFrame, de_down: pd.DataFrame, thre: float = 1.5, ignore_nonde: bool = False, obj: 'json object' = None):
        """ Evaluate a cluster to determine its putative cell type.
        """
        if obj is None:
            obj = self.object

        results = []
        for celltype in obj["cell_types"]:
            ct = CellType(celltype["name"], ignore_nonde = ignore_nonde)
            ct.evaluate(celltype, de_up, de_down, thre)
            if ct.score >= 0.5:
                sub_obj = celltype.get("subtypes", None)
                if sub_obj is not None:
                    ct.subtypes = self.evaluate(
                        de_up, de_down, thre = thre, ignore_nonde = ignore_nonde, obj = sub_obj
                    )
            results.append(ct)

        results.sort(key=lambda x: x.score, reverse=True)

        return results

    def report(self, fout: 'output stream', results, thre: float, space: int = 4) -> None:
        """ Write putative cell type reports to fout.
        """
        for ct in results:
            if ct.score >= thre:
                fout.write(" " * space + ct.tostring() + "\n")
                if ct.subtypes is not None:
                    self.report(fout, ct.subtypes, 0.5, space + 4)


def infer_cell_types(
    data: 'AnnData',
    markers: Union[str, Dict],
    de_test: str,
    de_alpha: float = 0.05,
    de_key: str = 'de_res',
    threshold: float = 0.5,
    ignore_nonde: bool = False,
    fout: 'output stream' = stdout,
) -> None:
    """Infer putative cell types for each cluster using legacy markers

    Parameters
    ----------

    data : `AnnData`
        AnnData object.
    markers : `str`
        A JSON file containing legacy markers or a python dictionary describing the markers. If you want to use predefined scCloud markers, use 'human_immune' for human immune cells, 'mouse_immune' for mouse immune cells, 'human_brain' for human brain cells, and 'mouse_brain' for mouse brain cells.
    de_test: `str`,
        scCloud determines cell types using DE test results. This argument indicates which DE test result to use, can be either 't', 'fisher' or 'mwu'.
    de_alpha: `float`, optional (default: 0.05)
        False discovery rate for controling family-wide error.
    de_key : `str`, optional (default: de_res)
        The keyword in varm that stores DE results.
    threshold : `float`, optional (defaut: 0.5)
        Only report putative cell types with a score larger than threshold.
    ignore_nonde: `bool`, optional (default: False)
        Do not consider non DE genes as weak negative markers.
    fout: `output stream`, optional (default: stdout)
        Where to write results. If default, write to stdout.

    Returns
    -------

    None

    Examples
    --------

    >>> annotate_cluster.infer_cell_types(adata, 'human_immune', 'fisher')
    """
    import pkg_resources
    if markers == "human_immune":
        markers = pkg_resources.resource_filename(
            "scCloud.annotate_cluster", "human_immune_cell_markers.json"
        )
    elif markers == "mouse_immune":
        markers = pkg_resources.resource_filename(
            "scCloud.annotate_cluster", "mouse_immune_cell_markers.json"
        )
    elif markers == "mouse_brain":
        markers = pkg_resources.resource_filename(
            "scCloud.annotate_cluster", "mouse_brain_cell_markers.json"
        )
    elif markers == "human_brain":
        markers = pkg_resources.resource_filename(
            "scCloud.annotate_cluster", "human_brain_cell_markers.json"
        )

    anno = Annotator(markers, data.var_names)

    clusts = natsorted([x.rpartition(':')[2] for x in data.varm[de_key].dtype.names if x.startswith("WAD_score:")])

    for clust_id in clusts:
        idx = data.varm[de_key]["{0}_qval:{1}".format(de_test, clust_id)] <= de_alpha

        idx_up = idx & (data.varm[de_key]["log_fold_change:{0}".format(clust_id)] > 0.0)
        idx_down = idx & (data.varm[de_key]["log_fold_change:{0}".format(clust_id)] < 0.0)
        assert idx_up.sum() + idx_down.sum() == idx.sum()

        cols = ["{0}:{1}".format(x, clust_id) for x in ["percentage_fold_change" if de_test == "fisher" else "log_fold_change", "percentage"]]
        de_up = pd.DataFrame(data = data.varm[de_key][cols][idx_up], index = data.var_names[idx_up])
        de_up.rename(columns={cols[0]: "fc", cols[1]: "percent"}, inplace=True)
        de_down = pd.DataFrame(data = data.varm[de_key][cols][idx_down], index = data.var_names[idx_down])
        de_down.rename(columns={cols[0]: "fc", cols[1]: "percent"}, inplace=True)

        if de_test != "fisher":
            de_up["fc"] = np.exp(de_up["fc"])
            de_down["fc"] = np.exp(de_down["fc"])

        results = anno.evaluate(de_up, de_down, ignore_nonde = ignore_nonde)

        fout.write("Cluster {0}:\n".format(clust_id))
        anno.report(fout, results, threshold)


def annotate(data: 'AnnData', name: str, based_on: str, anno_dict: Dict[str, str]) -> None:
    """Add annotation to AnnData obj.

    Parameters
    ----------

    data : `AnnData`
        AnnData object.
    name : `str`
        Name of the new annotation in data.obs.
    based_on : `str`
        Name of the attribute the cluster ids coming from.
    anno_dict : `Dict[str, str]`
        Dictionary mapping from cluster id to cell type.

    Returns
    -------

    None

    Examples
    --------
    >>> annotate_cluster.annotate(data, 'anno', 'spectral_louvain_labels', {'1': 'T cell', '2': 'B cell'})
    """
    data.obs[name] = [anno_dict[x] for x in data.obs[based_on]]


def run_annotate_cluster(
    input_file: str,
    output_file: str,
    marker_file: str,
    de_test: str,
    de_alpha: float = 0.05,
    de_key: str = 'de_res',
    threshold: float = 0.5,
    ignore_nonde: bool = False,
) -> None:
    """ For command line use.
    """
    import time
    from scCloud.io import read_input

    start = time.time()
    data = read_input(input_file, h5ad_mode="r")
    with open(output_file, "w") as fout:
        infer_cell_types(data, marker_file, de_test, de_alpha = de_alpha, de_key = de_key, threshold = threshold, ignore_nonde = ignore_nonde, fout = fout)
    data.file.close()
    end = time.time()
    logger.info("Time spent for annotating clusters is {:.2f}s.".format(end - start))


def annotate_anndata_object(input_file: str, annotation: str) -> None:
    """ For command line use.
        annotation:  anno_name:clust_name:cell_type1;...cell_typen
    """
    from scCloud.io import read_input, write_output

    data = read_input(input_file, h5ad_mode="r+")
    anno_name, clust_name, anno_str = annotation.split(":")
    anno_dict = {str(i + 1) : x for i, x in enumerate(anno_str.split(";"))}
    annotate(data, anno_name, clust_name, anno_dict)
    write_output(data, input_file)
