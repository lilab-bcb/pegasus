import numpy as np
import anndata

from typing import List
import logging

logger = logging.getLogger("pegasus")


def parse_subset_selections(subset_selections):
    subsets_dict = {}
    for subset_str in subset_selections:
        attr, value_str = subset_str.split(":")
        if attr in subsets_dict:
            subsets_dict[attr].extend(value_str.split(","))
        else:
            subsets_dict[attr] = value_str.split(",")
    return subsets_dict


def get_anndata_for_subclustering(data: "AnnData", subset_selections: List[str]):
    obs_index = np.full(data.shape[0], True)
    subsets_dict = parse_subset_selections(subset_selections)
    for key, value in subsets_dict.items():
        logger.info("{} in {}".format(str(key), str(value)))
        obs_index = obs_index & np.isin(data.obs[key], value)
    data = data[obs_index, :]

    obs_dict = {"obs_names": data.obs_names.values}
    for attr in data.obs.columns:
        if attr != "pseudotime":
            if attr.find("_labels") < 0:
                obs_dict[attr] = data.obs[attr].values
            else:
                obs_dict["parent_" + attr] = data.obs[attr].values

    var_dict = {
        "var_names": data.var_names.values,
        "gene_ids": data.var["gene_ids"].values,
        "robust": data.var["robust"].values,
    }

    newdata = anndata.AnnData(X=data.X, obs=obs_dict, var=var_dict)
    newdata.var["n_cells"] = newdata.X.getnnz(axis=0)
    newdata.var["robust"] = (
        newdata.var["robust"].values & (newdata.var["n_cells"] > 0).values
    )
    newdata.var["highly_variable_features"] = newdata.var[
        "robust"
    ]  # default all robust genes are "highly" variable

    if "Channels" in data.uns:
        newdata.uns["Channels"] = data.uns["Channels"]
    if "Groups" in data.uns:
        newdata.uns["Groups"] = data.uns["Groups"]
    if "plus" in data.varm.keys():
        newdata.varm["plus"] = data.varm["plus"]
    if "muls" in data.varm.keys():
        newdata.varm["muls"] = data.varm["muls"]

    logger.info("{0} cells are selected.".format(newdata.shape[0]))

    return newdata
