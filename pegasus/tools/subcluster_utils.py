import numpy as np
from pegasusio import MultimodalData
from typing import List, Dict, Union

import logging
logger = logging.getLogger(__name__)


def parse_subset_selections(subset_selections: List[str]) -> Dict[str, List[str]]:
    """ Take a list of selection criteria strings and return a parsed selection dictionary.

    Parameters
    ----------
    subset_selections: ``List[str]``
        A list of subset selection criteria. Each criterion takes the format of 'attribute:value,value,...,value'.

    Returns
    -------
    ``Dict[str, List[str]]``
    A dictionary of subsetting criteria
    """
    subsets_dict = {}
    for subset_str in subset_selections:
        attr, value_str = subset_str.split(":")
        if attr in subsets_dict:
            subsets_dict[attr].extend(value_str.split(","))
        else:
            subsets_dict[attr] = value_str.split(",")
    return subsets_dict


def clone_subset(data: MultimodalData, subset_selections: Union[List[bool], Union[str, List[str]]], ):
    """ Clone a subset of data as a new MultimodalData object for subclustering purpose.

    Parameters
    ----------
    subset_selections: ``List[bool]`` or ``str`` or ``List[str]``
        This parameter can take three forms. First, it can be an array/list of boolean values (List[bool]) indicating which cells are selected. Second, it can be a string with the format of 'attribute:value,value,...,value'. In this case, Pegasus will search for 'attribute' in the 'obs' field and select all cells with 'attribute' in the specified values. Lastly, it can be a list of strings as criteria and Pegasus will only select cells that satisfy all criteria (AND logic).

    Returns
    -------
    ``Dict[str, List[str]]``
    A dictionary of subsetting criteria
    """

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
