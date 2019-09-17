import numpy as np
import pandas as pd
import os
import time
from subprocess import check_call

from typing import List
from anndata import AnnData

from sccloud.io import infer_file_format, read_input, write_output, MemData


def find_digits(value):
    pos = len(value) - 1
    while pos >= 0 and value[pos].isdigit():
        pos -= 1
    pos += 1
    assert pos < len(value)
    return (value[:pos], int(value[pos:]))


def parse_restriction_string(rstr):
    pos = rstr.index(":")
    name = rstr[:pos]
    isin = True
    if rstr[pos + 1] == "~":
        isin = False
        pos += 1
    content = set()
    for item in rstr[pos + 1 :].split(","):
        values = item.split("-")
        if len(values) == 1:
            content.add(values[0])
        else:
            prefix, fr = find_digits(values[0])
            assert values[1].isdigit()
            to = int(values[1]) + 1
            for i in range(fr, to):
                content.add(prefix + str(i))
    return (name, isin, content)


def aggregate_matrices(
    csv_file: str,
    what_to_return: str = AnnData,
    restrictions: List[str] = [],
    attributes: List[str] = [],
    google_cloud: bool = False,
    select_singlets: bool = False,
    ngene: int = None,
    concat_matrices: bool = False,
) -> "None or AnnData or MemData":
    """Aggregate channel-specific count matrices into one big count matrix.

    This function takes as input a csv_file, which contains at least 2 columns — Sample, sample name; Location, file that contains the count matrices (e.g. filtered_gene_bc_matrices_h5.h5), and merges matrices from the same genome together. Depending on what_to_return, it can output the merged results into a sccloud-formatted HDF5 file or return as an AnnData or MemData object.

    Parameters
    ----------

    csv_file : `str`
        The CSV file containing information about each channel.
    what_to_return : `str`, optional (default: 'AnnData')
        If this value is equal to 'AnnData' or 'MemData', an AnnData or MemData object will be returned. Otherwise, results will be written into 'what_to_return.sccloud.h5' file and None is returned.
    restrictions : `list[str]`, optional (default: [])
        A list of restrictions used to select channels, each restriction takes the format of name:value,…,value or name:~value,..,value, where ~ refers to not.
    attributes : `list[str]`, optional (default: [])
        A list of attributes need to be incorporated into the output count matrix.
    google_cloud : `bool`, optional (default: False)
        If the channel-specific count matrices are stored in a google bucket.
    select_singlets : `bool`, optional (default: False)
        If we have demultiplexed data, turning on this option will make sccloud only include barcodes that are predicted as singlets.
    ngene : `int`, optional (default: None)
        The minimum number of expressed genes to keep one barcode.
    concat_matrices : `bool`, optional (default: False)
        If concatenate multiple matrices. If so, return only one AnnData object, otherwise, might return a list of AnnData objects.

    Returns
    -------

    None

    Examples
    --------
    >>> tools.aggregate_matrix('example.csv', 'example_10x.h5', ['Source:pbmc', 'Donor:1'], ['Source', 'Platform', 'Donor'])
    """

    df = pd.read_csv(csv_file, header=0, index_col="Sample")
    df["Sample"] = df.index

    # Select channels
    rvec = [parse_restriction_string(x) for x in restrictions]

    idx = pd.Series([True] * df.shape[0], index=df.index, name="selected")
    for name, isin, content in rvec:
        assert name in df.columns
        if isin:
            idx = idx & df[name].isin(content)
        else:
            idx = idx & (~(df[name].isin(content)))

    df = df.loc[idx]

    if df.shape[0] == 0:
        raise ValueError("No channels pass the restrictions!")

    # Load channels
    tot = 0
    aggrData = MemData()
    dest_paths = []
    for sample_name, row in df.iterrows():
        input_file = os.path.expanduser(
            os.path.expandvars(row["Location"].rstrip(os.sep))
        )
        file_format, copy_path, copy_type = infer_file_format(input_file)
        if google_cloud:
            base_name = os.path.basename(copy_path)
            dest_path = sample_name + "_tmp_" + base_name

            if copy_type == "directory":
                check_call(["mkdir", "-p", dest_path])
                call_args = ["gsutil", "-m", "cp", "-r", copy_path, dest_path]
            else:
                call_args = ["gsutil", "-m", "cp", copy_path, dest_path]
            check_call(call_args)
            dest_paths.append(dest_path)

            input_file = dest_path
            if file_format == "csv" and copy_type == "directory":
                input_file = os.path.join(dest_path, os.path.basename(input_file))

        genome = None
        if file_format in ["dge", "csv", "mtx", "loom"]:
            assert "Reference" in row
            genome = row["Reference"]

        data = read_input(
            input_file,
            genome=genome,
            return_type="MemData",
            ngene=ngene,
            select_singlets=select_singlets,
        )
        data.update_barcode_metadata_info(sample_name, row, attributes)
        aggrData.addAggrData(data)

        tot += 1
        print("Processed {}.".format(input_file))

    # Delete temporary file
    for dest_path in dest_paths:
        check_call(["rm", "-rf", dest_path])

    # Merge channels
    t1 = time.time()
    aggrData.aggregate()
    t2 = time.time()
    print("Data aggregation is finished in {:.2f}s.".format(t2 - t1))

    if what_to_return == "AnnData":
        aggrData = aggrData.convert_to_anndata(concat_matrices)
    elif what_to_return != "MemData":
        write_output(aggrData, what_to_return)
        aggrData = None

    print("Aggregated {tot} files.".format(tot=tot))

    return aggrData
