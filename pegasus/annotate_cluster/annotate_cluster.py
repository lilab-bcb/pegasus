import numpy as np
import pandas as pd
import json

from sys import stdout
from natsort import natsorted
from typing import List, Dict, Union
from anndata import AnnData
from io import IOBase

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
        obj: dict,
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
                                (marker, f"{percent:.2f}%")
                            )
                        else:
                            numer += 1.0 + (fc - 1.0) / (thre - 1.0)
                            self.weak_support.append(
                                (marker, f"{percent:.2f}%")
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
                                    (marker, f"{percent:.2f}%")
                                )
                            else:
                                numer += 1.0 + (fc - 1.0) / (thre - 1.0)
                                self.weak_support.append(
                                    (marker, f"{percent:.2f}%")
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
        res = f"name: {self.name}; score: {self.score:.2f}; average marker percentage: {self.avgp:.2f}%"
        if len(self.strong_support) > 0:
            res += "; strong support: {0}".format(
                ",".join([f"({x[0]},{x[1]})" for x in self.strong_support])
            )
        if len(self.weak_support) > 0:
            res += "; weak support: {0}".format(
                ",".join([f"({x[0]},{x[1]})" for x in self.weak_support])
            )

        return res


class Annotator:
    def __init__(self, markers: Dict, genes: List[str]) -> None:
        self.object = markers
        self.recalibrate(self.object, genes)

    def recalibrate(self, obj: dict, genes: List[str]) -> None:
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
        fc_thre: float = 1.5,
        threshold: float = 0.5,
        ignore_nonde: bool = False,
        obj: dict = None,
    ):
        """ Evaluate a cluster to determine its putative cell type.
        """
        if obj is None:
            obj = self.object

        results = []
        for celltype in obj["cell_types"]:
            ct = CellType(celltype["name"], ignore_nonde=ignore_nonde)
            ct.evaluate(celltype, de_up, de_down, fc_thre)
            if ct.score >= threshold:
                sub_obj = celltype.get("subtypes", None)
                if sub_obj is not None:
                    ct.subtypes = self.evaluate(
                        de_up,
                        de_down,
                        fc_thre=fc_thre,
                        ignore_nonde=ignore_nonde,
                        obj=sub_obj,
                    )
                results.append(ct)

        results.sort(key=lambda x: x.score, reverse=True)

        return results

    def report(
        self,
        fout: IOBase,
        ct_list: List["CellType"],
        space: int = 4,
    ) -> None:
        """ Write putative cell type reports to fout.
        """
        for ct in ct_list:
            fout.write(" " * space + str(ct) + "\n")
            if ct.subtypes is not None:
                self.report(fout, ct.subtypes, space + 4)


def infer_cluster_names(
    cell_type_dict: Dict[str, List["CellType"]], threshold: float = 0.5, is_human_immune: bool = False
) -> List[str]:
    """Decide cluster names based on cell types automatically.

    Parameters
    ----------
    cell_type_dict: ``Dict[str, List["CellType"]]``
        Python dictionary of cell type list for each cluster. This is the output of ``pg.infer_cell_types``.

    threshold: ``float``, optional, default: ``0.5``
        A threshold for cell type result reported. It should be a real number between ``0.0`` and ``1.0``.

    is_human_immune: ``bool``, optional, default: False
        If cell types are annotated using the Pegasus' human_immune markers, apply a smarter naming algorithm.

    Returns
    -------
    ``List[str]``
        A list of cluster names decided by their corresponding cell types. The order is consistent with clusters.

    Examples
    --------
    >>> cell_type_dict = pg.infer_cell_types(adata, markers = 'human_immune', de_test = 't')
    >>> cluster_names = pg.infer_cluster_names(cell_type_dict)
    """
    from collections import Counter

    cluster_ids = natsorted(cell_type_dict.keys())
    names = []
    name_dict = Counter()
    for cluster_id in cluster_ids:
        ct_list = cell_type_dict[cluster_id]

        if len(ct_list) == 0 or ct_list[0].score < threshold:
            cell_name = cluster_id
        else:
            ct = ct_list[0]
            if is_human_immune and ct.name == "T cell":
                subname = None
                has_naive_t = False
                for subt in ct.subtypes:
                    if subt.score >= threshold and (subt.name != "T regulatory cell" or subt.avgp > 0.5):
                        if subt.name == "Naive T cell" and subt.score >= 0.6:
                            has_naive_t = True
                        elif subname is None:
                            subname = subt.name
                if subname is None:
                    cell_name = "Naive T cell" if has_naive_t else "T cell"
                elif has_naive_t and (subname in ["T helper cell", "Cytotoxic T cell"]):
                    cell_name = "CD4+ Naive T cell" if subname == "T helper cell" else "CD8+ Naive T cell"
                else:
                    cell_name = subname
            elif is_human_immune and ct.name == "B cell":
                subname = None
                for subt in ct.subtypes:
                    if subt.score >= threshold:
                        subname = subt.name
                        break
                cell_name = subname if subname is not None else "Naive B cell"
            elif is_human_immune and ct.name == "CD1C+ dendritic cell":
                cell_name = ct.name
                for ctype in ct_list[1:]:
                    if ctype.score >= threshold and ctype.name == "CLEC9A+ dendritic cell":
                        cell_name = "Conventional dendritic cell (CD1C+/CLEC9A+)"
                        break
            else:
                while ct.subtypes is not None and len(ct.subtypes) > 0 and ct.subtypes[0].score >= threshold:
                    ct = ct.subtypes[0]
                cell_name = ct.name

            name_dict[cell_name] += 1
            if name_dict[cell_name] > 1:
                cell_name = f"{cell_name}-{name_dict[cell_name]}"

        names.append(cell_name)

    return names


def infer_cell_types(
    data: Union[MultimodalData, UnimodalData, AnnData],
    markers: Union[str, Dict],
    de_test: str = "mwu",
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
        * If ``str``, it is a string representing a comma-separated list; each element in the list
            * either refers to a JSON file containing legacy markers, or predefined markers
            * ``'human_immune'`` for human immune cells;
            * ``'mouse_immune'`` for mouse immune cells;
            * ``'human_brain'`` for human brain cells;
            * ``'mouse_brain'`` for mouse brain cells;
            * ``'human_lung'`` for human lung cells.
        * If ``Dict``, it refers to a Python dictionary describing the markers.

    de_test: ``str``, optional, default: ``"mwu"``
        pegasus determines cell types using DE test results. This argument indicates which DE test result to use, can be either ``'t'``, ``'fisher'`` or ``'mwu'``.
        By default, it uses ``'mwu'``.

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
    >>> cell_type_dict = pg.infer_cell_types(adata, markers = 'human_immune,human_brain')
    """

    if output_file is not None:
        fout = open(output_file, "w")

    import pkg_resources
    predefined_markers = dict(
        human_immune="human_immune_cell_markers.json",
        mouse_immune="mouse_immune_cell_markers.json",
        human_brain="human_brain_cell_markers.json",
        mouse_brain="mouse_brain_cell_markers.json",
        human_lung="human_lung_cell_markers.json",
    )

    if isinstance(markers, str):
        tokens = markers.split(',')
        markers = None
        for token in tokens:
            if token in predefined_markers:
                token = pkg_resources.resource_filename(
                    "pegasus.annotate_cluster", predefined_markers[token]
                )
            with open(token) as fin:
                tmp_dict = json.load(fin)
            if markers is None:
                markers = tmp_dict
            else:
                markers["title"] = f"{markers['title']}/{tmp_dict['title']}"
                markers["cell_types"].extend(tmp_dict["cell_types"])

    assert isinstance(markers, dict)
    anno = Annotator(markers, data.var_names)

    test2metric = {"mwu": "auroc", "t": "log2FC", "fisher": "percentage_fold_change"}
    metric = test2metric[de_test]
    thre = 0.5 if de_test == "mwu" else 0.0
    coln = "percentage_fold_change" if de_test == "fisher" else "log2FC"

    clusts = natsorted(
        [
            x.rpartition(":")[0]
            for x in data.varm[de_key].dtype.names
            if x.endswith(":auroc")
        ]
    )
    cell_type_results = {}
    for clust_id in clusts:
        idx = data.varm[de_key][f"{clust_id}:{de_test}_qval"] <= de_alpha

        idx_up = idx & (data.varm[de_key][f"{clust_id}:{metric}"] > thre)
        idx_down = idx & (data.varm[de_key][f"{clust_id}:{metric}"] < thre)
        assert idx_up.sum() + idx_down.sum() == idx.sum()

        cols = [f"{clust_id}:{coln}", f"{clust_id}:percentage"]
        de_up = pd.DataFrame(
            data=data.varm[de_key][cols][idx_up], index=data.var_names[idx_up]
        )
        de_up.rename(columns={cols[0]: "fc", cols[1]: "percent"}, inplace=True)
        de_down = pd.DataFrame(
            data=data.varm[de_key][cols][idx_down], index=data.var_names[idx_down]
        )
        de_down.rename(columns={cols[0]: "fc", cols[1]: "percent"}, inplace=True)

        if de_test != "fisher":
            de_up["fc"] = 2.0 ** de_up["fc"]
            de_down["fc"] = 2.0 ** de_down["fc"]

        results = anno.evaluate(de_up, de_down, threshold=threshold, ignore_nonde=ignore_nonde)

        if output_file is not None:
            fout.write(f"Cluster {clust_id}:\n")
            anno.report(fout, results)

        cell_type_results[clust_id] = results

    if output_file is not None:
        fout.close()

    return cell_type_results


def annotate(
    data: Union[MultimodalData, UnimodalData,AnnData],
    name: str,
    based_on: str,
    anno_dict: Union[Dict[str, str], List[str]],
) -> None:
    """Add annotation to the data object as a categorical variable.

    Parameters
    ----------

    data : ``MultimodalData``, ``UnimodalData``, or ``anndata.AnnData``
        Gene-count matrix with DE analysis information.
    name : `str`
        Name of the new annotation in data.obs.
    based_on : `str`
        Name of the attribute the cluster ids coming from.
    anno_dict : `Dict[str, str]` or `List[str]`
        Dictionary mapping from cluster id to cell type.
        If it is a List, map cell types to cluster ids one to one in correspondence.

    Returns
    -------

    None

    Examples
    --------
    >>> pg.annotate(data, 'anno', 'spectral_louvain_labels', {'1': 'T cell', '2': 'B cell'})
    >>> pg.annotate(data, 'anno', 'louvain_labels', ['T cell', 'B cell'])
    """
    if isinstance(anno_dict, list):
        cluster_ids = data.obs[based_on].cat.categories.values.astype('str')
        anno_dict = dict(zip(cluster_ids, anno_dict))
    from natsort import natsorted 
    data.obs[name] = pd.Categorical([anno_dict[x] for x in data.obs[based_on]], categories = natsorted(np.unique(list(anno_dict.values()))))

@timer(logger=logger)
def run_annotate_cluster(
    input_file: str,
    output_file: str,
    markers: str,
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
        markers,
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
