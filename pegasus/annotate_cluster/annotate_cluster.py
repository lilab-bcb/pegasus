import numpy as np
import pandas as pd
import json

from sys import stdout
from natsort import natsorted
from typing import List, Dict, Union
from anndata import AnnData

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer, MultimodalData, UnimodalData


class CellType:
    def __init__(self, name: str, ignore_nonde: bool = False):
        self.name = name
        self.score = self.avgp = 0.0
        self.weak_support = []
        self.strong_support = []
        self.subtypes = None
        self.ignore_nonde = ignore_nonde

    def evaluate(
        self,
        obj: "json object",
        de_up: pd.DataFrame,
        de_down: pd.DataFrame,
        thre: float,
    ):
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
                            self.strong_support.append(
                                (marker, "{0:.2f}%".format(percent))
                            )
                        else:
                            numer += 1.0 + (fc - 1.0) / (thre - 1.0)
                            self.weak_support.append(
                                (marker, "{0:.2f}%".format(percent))
                            )
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

    def __repr__(self):
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

    def recalibrate(self, obj: "json object", genes: List[str]) -> None:
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

    def evaluate(
        self,
        de_up: pd.DataFrame,
        de_down: pd.DataFrame,
        thre: float = 1.5,
        ignore_nonde: bool = False,
        obj: "json object" = None,
    ):
        """ Evaluate a cluster to determine its putative cell type.
        """
        if obj is None:
            obj = self.object

        results = []
        for celltype in obj["cell_types"]:
            ct = CellType(celltype["name"], ignore_nonde=ignore_nonde)
            ct.evaluate(celltype, de_up, de_down, thre)
            if ct.score >= 0.5:
                sub_obj = celltype.get("subtypes", None)
                if sub_obj is not None:
                    ct.subtypes = self.evaluate(
                        de_up,
                        de_down,
                        thre=thre,
                        ignore_nonde=ignore_nonde,
                        obj=sub_obj,
                    )
            results.append(ct)

        results.sort(key=lambda x: x.score, reverse=True)

        return results

    def report(
        self,
        fout: "output stream",
        ct_list: List["CellType"],
        thre: float,
        space: int = 4,
    ) -> None:
        """ Write putative cell type reports to fout.
        """
        for ct in ct_list:
            if ct.score >= thre:
                fout.write(" " * space + str(ct) + "\n")
                if ct.subtypes is not None:
                    self.report(fout, ct.subtypes, 0.5, space + 4)


def infer_cluster_names(
    cell_type_dict: Dict[str, List["CellType"]], threshold: float = 0.5
) -> List[str]:
    """Decide cluster names based on cell types automatically.

    Parameters
    ----------
    cell_type_dict: ``Dict[str, List["CellType"]]``
        Python dictionary of cell type list for each cluster. This is the output of ``pg.infer_cell_types``.

    threshold: ``float``, optional, default: ``0.5``
        A threshold for cell type result reported. It should be a real number between ``0.0`` and ``1.0``.

    Returns
    -------
    ``List[str]``
        A list of cluster names decided by their corresponding cell types. The order is consistent with clusters.

    Examples
    --------
    >>> cell_type_dict = pg.infer_cell_types(adata, markers = 'human_immune', de_test = 't')
    >>> cluster_names = pg.infer_cluster_names(cell_type_dict)
    """
    cluster_ids = natsorted(cell_type_dict.keys())
    names = []
    name_dict = dict()
    for cluster_id in cluster_ids:
        ct_list = cell_type_dict[cluster_id]

        if len(ct_list) == 0 or ct_list[0].score < threshold:
            cell_name = cluster_id
        else:
            ct = ct_list[0]
            while ct.subtypes is not None and ct.subtypes[0].score >= threshold:
                ct = ct.subtypes[0]
            cell_name = ct.name

            if cell_name in name_dict:
                name_dict[cell_name] += 1
                cell_name = cell_name + "-" + str(name_dict[cell_name])
            else:
                name_dict[cell_name] = 1

        names.append(cell_name)

    return names


def infer_cell_types(
    data: Union[MultimodalData, UnimodalData, AnnData],
    markers: Union[str, Dict],
    de_test: str,
    de_alpha: float = 0.05,
    de_key: str = "de_res",
    threshold: float = 0.5,
    ignore_nonde: bool = False,
    output_file: str = None,
) -> Dict[str, List["CellType"]]:
    """Infer putative cell types for each cluster using legacy markers.

    Parameters
    ----------

    data : ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``.
        Data structure of count matrix and DE analysis results.

    markers : ``str`` or ``Dict``
        * If ``str``, it
            * either refers to a JSON file containing legacy markers, or
            * ``'human_immune'`` for predefined pegasus markers on human immune cells;
            * ``'mouse_immune'`` for mouse immune cells;
            * ``'human_brain'`` for human brain cells;
            * ``'mouse_brain'`` for mouse brain cells.
        * If ``Dict``, it refers to a Python dictionary describing the markers.

    de_test: ``str``
        pegasus determines cell types using DE test results. This argument indicates which DE test result to use, can be either ``'t'``, ``'fisher'`` or ``'mwu'``.

    de_alpha: ``float``, optional, default: ``0.05``
        False discovery rate for controling family-wide error.

    de_key : ``str``, optional, default: ``"de_res"``
        The keyword in ``data.varm`` that stores DE analysis results.

    threshold : ``float``, optional, defaut: ``0.5``
        Only report putative cell types with a score larger than or equal to ``threshold``.

    ignore_nonde: ``bool``, optional, default: ``False``
        Do not consider non DE genes as weak negative markers.

    output_file: ``str``, optional, default: ``None``
        File name of output cluster annotation. If ``None``, do not write to any file.

    Returns
    -------
    ``Dict[str, List["CellType"]]``
        Python dictionary with cluster ID's being keys, and their corresponding cell type lists sortec by scores being values.

    Examples
    --------
    >>> cell_type_dict = pg.infer_cell_types(adata, markers = 'human_immune', de_test = 't')
    """

    if output_file is not None:
        fout = open(output_file, "w")

    predefined_markers = dict(
        human_immune="human_immune_cell_markers.json",
        mouse_immune="mouse_immune_cell_markers.json",
        human_brain="human_brain_cell_markers.json",
        mouse_brain="mouse_brain_cell_markers.json",
        human_lung="human_lung_cell_markers.json",
    )

    if markers in predefined_markers:
        import pkg_resources

        markers = pkg_resources.resource_filename(
            "pegasus.annotate_cluster", predefined_markers[markers]
        )

    anno = Annotator(markers, data.var_names)

    clusts = natsorted(
        [
            x.rpartition(":")[2]
            for x in data.varm[de_key].dtype.names
            if x.startswith("WAD_score:")
        ]
    )
    cell_type_results = {}
    for clust_id in clusts:
        idx = data.varm[de_key]["{0}_qval:{1}".format(de_test, clust_id)] <= de_alpha

        idx_up = idx & (data.varm[de_key]["log_fold_change:{0}".format(clust_id)] > 0.0)
        idx_down = idx & (
            data.varm[de_key]["log_fold_change:{0}".format(clust_id)] < 0.0
        )
        assert idx_up.sum() + idx_down.sum() == idx.sum()

        cols = [
            "{0}:{1}".format(x, clust_id)
            for x in [
                "percentage_fold_change" if de_test == "fisher" else "log_fold_change",
                "percentage",
            ]
        ]
        de_up = pd.DataFrame(
            data=data.varm[de_key][cols][idx_up], index=data.var_names[idx_up]
        )
        de_up.rename(columns={cols[0]: "fc", cols[1]: "percent"}, inplace=True)
        de_down = pd.DataFrame(
            data=data.varm[de_key][cols][idx_down], index=data.var_names[idx_down]
        )
        de_down.rename(columns={cols[0]: "fc", cols[1]: "percent"}, inplace=True)

        if de_test != "fisher":
            de_up["fc"] = np.exp(de_up["fc"])
            de_down["fc"] = np.exp(de_down["fc"])

        results = anno.evaluate(de_up, de_down, ignore_nonde=ignore_nonde)

        if output_file is not None:
            fout.write("Cluster {}:\n".format(clust_id))
            anno.report(fout, results, threshold)

        cell_type_results[clust_id] = results

    if output_file is not None:
        fout.close()

    return cell_type_results


def annotate(
    data: Union[MultimodalData, UnimodalData,AnnData],
    name: str,
    based_on: str,
    anno_dict: Dict[str, str],
) -> None:
    """Add annotation to AnnData obj.

    Parameters
    ----------

    data : ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``
        Gene-count matrix with DE analysis information.
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

@timer(logger=logger)
def run_annotate_cluster(
    input_file: str,
    output_file: str,
    marker_file: str,
    de_test: str,
    de_alpha: float = 0.05,
    de_key: str = "de_res",
    threshold: float = 0.5,
    ignore_nonde: bool = False,
) -> None:
    """ For command line use.
    """
    from pegasusio import read_input

    data = read_input(input_file, mode="r")
    infer_cell_types(
        data,
        marker_file,
        de_test,
        de_alpha=de_alpha,
        de_key=de_key,
        threshold=threshold,
        ignore_nonde=ignore_nonde,
        output_file=output_file,
    )


def annotate_data_object(input_file: str, annotation: str) -> None:
    """ For command line use.
        annotation:  anno_name:clust_name:cell_type1;...cell_typen
    """
    from pegasusio import read_input, write_output

    data = read_input(input_file, mode="r")
    anno_name, clust_name, anno_str = annotation.split(":")
    anno_dict = {str(i + 1): x for i, x in enumerate(anno_str.split(";"))}
    annotate(data, anno_name, clust_name, anno_dict)
    write_output(data, input_file)
