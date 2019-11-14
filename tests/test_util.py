import os
import subprocess

import numpy as np
import pandas as pd
import pegasus as pg
import scipy


def is_running_in_docker():
    return os.path.exist("/.dockerenv")


def assert_excel_equal(test_case, path, test_path):
    test_output = pd.read_excel(test_path, sheet_name=None)
    output = pd.read_excel(path, sheet_name=None)
    errors = []
    try:
        test_case.assertListEqual(list(test_output.keys()), list(output.keys()))
    except AssertionError as ex:
        errors.append(str(ex))
    for key in output:
        errors += assert_data_frames_equal(test_case, test_output[key], output[key])
    test_case.assertEqual(len(errors), 0, "\n".join(errors))


def assert_files_equal(test_case, path1, path2):
    c = subprocess.run(["diff", "-q", path1, path2])
    test_case.assertEqual(c.returncode, 0, "{} and {} differ".format(path1, path2))


def assert_dict_of_arrays_equal(test_case, dict1, dict2, blacklist):
    # TODO handle nested
    if blacklist is None:
        blacklist = set()
    keys1 = set(dict1.keys())
    keys2 = set(dict2.keys())
    errors = []
    try:
        test_case.assertSetEqual(
            keys1.difference(blacklist), keys2.difference(blacklist)
        )
    except AssertionError as ex:
        errors.append(str(ex))
    all_keys = set(list(dict1.keys()) + list(dict2.keys()))
    for key in all_keys:
        if key not in blacklist and key in keys1 and key in keys2:
            val1 = dict1[key]
            val2 = dict2[key]
            if scipy.sparse.issparse(val1):
                val1 = val1.toarray()
            if scipy.sparse.issparse(val2):
                val2 = val2.toarray()
            try:
                if isinstance(val1, np.ndarray):
                    np.testing.assert_array_almost_equal(
                        val1[~np.isnan(val1)],
                        val2[~np.isnan(val2)],
                        err_msg=str(key) + " not equal",
                    )
                else:
                    test_case.assertEqual(val1, val2, str(key) + " not equal")
                    # np.testing.assert_array_equal(val1[~np.isnan(val1)], val2[~np.isnan(val2)],
                    #     err_msg=str(key) + ' not equal')
            except AssertionError as ex:
                errors.append(str(ex))
    return errors


def assert_data_frames_equal(test_case, df1, df2, blacklist=None):
    if blacklist is None:
        blacklist = set()
    list1_keys = set(df1.columns)
    list2_keys = set(df2.columns)
    errors = []
    try:
        test_case.assertSetEqual(
            list1_keys.difference(blacklist), list2_keys.difference(blacklist)
        )
    except AssertionError as ex:
        errors.append(str(ex))
    all_keys = set(list(df1.columns) + list(df2.columns))
    for key in all_keys:
        if key not in blacklist and key in df1 and key in df2:
            val1 = df1[~df1[key].isna()][key].values
            val2 = df2[~df2[key].isna()][key].values

            try:
                if val1.dtype.kind != "O":
                    np.testing.assert_array_almost_equal(
                        val1[~np.isnan(val1)],
                        val2[~np.isnan(val2)],
                        err_msg=str(key) + " not equal",
                    )
                else:
                    np.testing.assert_array_equal(
                        val1, val2, err_msg=str(key) + " not equal"
                    )
            except AssertionError as ex:
                errors.append(str(ex))
    return errors


def assert_adata_files_equal(
    test_case,
    path,
    test_path,
    obs_blacklist=None,
    var_blacklist=None,
    uns_blacklist=None,
    obsm_blacklist=None,
    varm_blacklist=None,
):
    test_data = pg.io.read_input(test_path)
    data = pg.io.read_input(path)
    assert_adata_equal(
        test_case,
        data,
        test_data,
        obs_blacklist,
        var_blacklist,
        uns_blacklist,
        obsm_blacklist,
        varm_blacklist,
    )


def assert_adata_equal(
    test_case,
    data1,
    data2,
    obs_blacklist=None,
    var_blacklist=None,
    uns_blacklist=None,
    obsm_blacklist=None,
    varm_blacklist=None,
):
    if data1.isbacked:
        data1.X = data1.X[()]
    if data2.isbacked:
        data2.X = data2.X[()]

    if scipy.sparse.issparse(data1.X):
        data1.X = data1.X.toarray()
    if scipy.sparse.issparse(data2.X):
        data2.X = data2.X.toarray()
    errors = []
    try:
        np.testing.assert_array_equal(data1.X, data2.X, err_msg="X not equal")
    except AssertionError as ex:
        errors.append(str(ex))
    if obs_blacklist is not None:
        obs1 = data1.obs.drop(obs_blacklist, axis=1)
        obs2 = data2.obs.drop(obs_blacklist, axis=1)
    else:
        obs1 = data1.obs
        obs2 = data2.obs

    errors += assert_data_frames_equal(test_case, obs1, obs2)
    if var_blacklist is not None:
        var1 = data1.var.drop(var_blacklist, axis=1)
        var2 = data2.var.drop(var_blacklist, axis=1)
    else:
        var1 = data1.var
        var2 = data2.var
    errors += assert_data_frames_equal(test_case, var1, var2)
    errors += assert_dict_of_arrays_equal(
        test_case, data1.uns, data2.uns, uns_blacklist
    )
    errors += assert_dict_of_arrays_equal(
        test_case, data1.obsm, data2.obsm, obsm_blacklist
    )
    errors += assert_dict_of_arrays_equal(
        test_case, data1.varm, data2.varm, varm_blacklist
    )
    test_case.assertEqual(len(errors), 0, "\n".join(errors))
