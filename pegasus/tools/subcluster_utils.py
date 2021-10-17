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


def clone_subset(data: MultimodalData, subset_selections: Union[List[bool], Union[str, List[str]]]) -> MultimodalData:
    """ Clone a subset of data as a new MultimodalData object for subclustering purpose.

    Parameters
    ----------
    data: ``pegasusio.MultimodalData``
        Annotated data matrix with rows for cells and columns for genes.
    subset_selections: ``List[bool]`` or ``str`` or ``List[str]``
        This parameter can take three forms. First, it can be an array/list of boolean values (List[bool]) indicating which cells are selected. Second, it can be a string with the format of 'attribute:value,value,...,value'. In this case, Pegasus will search for 'attribute' in the 'obs' field and select all cells with 'attribute' in the specified values. Lastly, it can be a list of strings as criteria and Pegasus will only select cells that satisfy all criteria (AND logic).

    Returns
    -------
    ``pegasusio.MultimodalData``
        Return a new subsetted data
    """
    if isinstance(subset_selections, str):
        subset_selections = [subset_selections]

    if len(subset_selections) == 0:
        raise ValueError(f"No entry in subset_selections!")

    val = subset_selections[0]
    if isinstance(val, bool) or isinstance(val, np.bool_): # if bool list
        obs_index = subset_selections
    else:
        obs_index = np.full(data.shape[0], True)
        subsets_dict = parse_subset_selections(subset_selections)
        for key, value in subsets_dict.items():
            logger.info(f"Select {key} in {str(value)}")
            obs_index = obs_index & np.isin(data.obs[key], value)

    newdata = data[obs_index, :].copy()
    # Add _par suffix to clustering labels
    mapper = {}
    for attr in newdata.obs.columns:
        if newdata.get_attr_type(attr) == "cluster":
            new_attr = f"{attr}_par"
            newdata.register_attr(attr)
            newdata.register_attr(new_attr, "cluster")
            mapper[attr] = new_attr
    if len(mapper) > 0:
        newdata.obs.rename(columns = mapper, inplace = True)
    # Add _par suffix to basis and also drop knn 
    for attr in list(newdata.obsm.keys()):
        attr_type = newdata.get_attr_type(attr)
        if attr_type == "basis":
            new_attr = f"{attr}_par"
            newdata.register_attr(attr)
            newdata.register_attr(new_attr, "basis")
            newdata.obsm[new_attr] = newdata.obsm.pop(attr)
        elif attr_type == "knn":
            newdata.register_attr(attr) # pop up
            del newdata.obsm[attr]        
    # For var, update 'n_cells', mark genes with n_cells == 0 as unrobust and make hvg same as robust by default
    newdata.var["n_cells"] = newdata.X.getnnz(axis=0)
    values = newdata.var["robust"].values if ("robust" in newdata.var) else np.full(newdata.shape[1], True)
    values &= (newdata.var["n_cells"] > 0).values
    newdata.var["robust"] = values
    newdata.var["highly_variable_features"] = newdata.var["robust"]  # default all robust genes are "highly" variable

    logger.info(f"{newdata.shape[0]} cells are selected.")

    return newdata
