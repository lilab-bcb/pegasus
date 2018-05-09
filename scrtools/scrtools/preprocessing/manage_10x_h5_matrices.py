#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from scipy.sparse import csc_matrix, hstack
import tables
import copy
from subprocess import check_call


def load_10x_h5_file(input_h5, genome):
    with tables.open_file(input_h5) as h5_in:
        attrs = {}
        for node in h5_in.walk_nodes("/" + genome, "Array"):
            attrs[node.name] = node.read()

        attrs["matrix"] = csc_matrix((attrs["data"], attrs["indices"], attrs["indptr"]), shape=attrs["shape"])

        attrs.pop("data")
        attrs.pop("indices")
        attrs.pop("indptr")
        attrs.pop("shape")

    return attrs


def write_10x_h5_file(output_h5, output_data, genome):
    with tables.open_file(output_h5, mode="w", title=output_h5, filters=tables.Filters(complevel=1)) as hd5_out:
        out_group = hd5_out.create_group("/", genome)
        hd5_out.create_carray(out_group, "barcodes", obj=output_data["barcodes"])
        hd5_out.create_carray(out_group, "gene_names", obj=output_data["gene_names"])
        hd5_out.create_carray(out_group, "genes", obj=output_data["genes"])
        hd5_out.create_carray(out_group, "data", obj=output_data["matrix"].data)
        hd5_out.create_carray(out_group, "indices", obj=output_data["matrix"].indices)
        hd5_out.create_carray(out_group, "indptr", obj=output_data["matrix"].indptr)
        hd5_out.create_carray(out_group, "shape", obj=output_data["matrix"].shape)


def find_digits(value):
    pos = len(value) - 1
    while pos >= 0 and value[pos].isdigit():
        pos -= 1
    pos += 1
    assert pos < len(value)
    return (value[:pos], int(value[pos:]))


def parse_restriction_string(rstr):
    pos = rstr.index(':')
    name = rstr[: pos]
    content = set()
    for item in rstr[pos + 1:].split(','):
        values = item.split('-')
        if len(values) == 1:
            content.add(values[0])
        else:
            prefix, fr = find_digits(values[0])
            assert values[1].isdigit()
            to = int(values[1]) + 1
            for i in range(fr, to):
                content.add(prefix + str(i))
    return (name, content)


def aggregate_10x_matrices(csv_file, genome, restrictions, attributes, output_name, google_cloud=False):
    df = pd.read_csv(csv_file, header=0, index_col='Sample')
    restrictions.append("Reference:{}".format(genome))
    rvec = [parse_restriction_string(x) for x in restrictions]
    idx = pd.Series([True] * df.shape[0], index=df.index, name='selected')
    for name, content in rvec:
        assert name in df.columns
        idx = idx & df[name].isin(content)
    df = df.loc[idx]

    if df.shape[0] == 0:
        print("No channels pass the restrictions!")
        return

    df[attributes].to_csv("{output_name}.attr.csv".format(output_name=output_name))

    tot = 0
    out_hd5 = None
    for sample_name, row in df.iterrows():
        input_hd5_file = "{location}/{sample}/filtered_gene_bc_matrices_h5.h5".format(location=row['Location'],
                                                                                      sample=sample_name)
        if google_cloud:
            call_args = ['gsutil', '-m', 'cp', input_hd5_file,
                         '{0}_filtered_gene_bc_matrices_h5.h5'.format(sample_name)]
            check_call(call_args)
            input_hd5_file = '{0}_filtered_gene_bc_matrices_h5.h5'.format(sample_name)
        elif not os.path.exists(input_hd5_file):
            input_hd5_file = "{location}/filtered_gene_bc_matrices_h5.h5".format(location=row['Location'])

        channel = load_10x_h5_file(input_hd5_file, genome)
        channel["barcodes"] = [sample_name + '-' + x.decode() for x in channel["barcodes"]]

        if out_hd5 is None:
            out_hd5 = copy.copy(channel)
            out_hd5["matrix"] = []
            out_hd5["barcodes"] = []

        out_hd5["matrix"].append(channel["matrix"])
        out_hd5["barcodes"].append(channel["barcodes"])

        print("Processed {}.".format(input_hd5_file))
        tot += 1

    out_hd5["matrix"] = hstack(out_hd5["matrix"], "csc")
    out_hd5["barcodes"] = [barcode for channel in out_hd5["barcodes"] for barcode in channel]

    write_10x_h5_file("{output_name}_10x.h5".format(output_name=output_name), out_hd5, genome)
    print("Generated {name}_10x.h5 and {name}.attr.csv from {tot} files.".format(name=output_name, tot=tot))


def add_attribute(attr_file, attribute_string):
    df = pd.read_csv(attr_file, header=0, index_col=0, dtype=str)
    new_attr, stri = attribute_string.split(':')
    assert new_attr not in df.columns
    if stri.find('=') >= 0:
        attr, value = stri.split('=')
        assert attr in df.columns
        df[new_attr] = (df[attr] == value).astype(int)
    elif stri.find('+') >= 0:
        attrs = stri.split('+')
        assert np.isin(attrs, df.columns).sum() == len(attrs)
        df[new_attr] = df[attrs].apply(lambda x: '+'.join(x), axis=1)
    else:
        df[new_attr] = stri
    df.to_csv(attr_file)


def merge_10x_matrices(input_names, rep_syms, genome, attributes, output_name):
    """ We append rep_syms in attributes """
    assert len(attributes) == 0 or len(rep_syms) == len(input_names)

    out_df = []
    out_hd5 = None
    for i, input_name in enumerate(input_names):
        input_hd5_file = "{input_name}_10x.h5".format(input_name=input_name)
        channel = load_10x_h5_file(input_hd5_file, genome)

        input_csv_file = "{input_name}.attr.csv".format(input_name=input_name)
        df = pd.read_csv(input_csv_file, header=0, index_col=0, dtype=str)
        for attr in attributes:
            assert attr in df.columns
            df[attr] = df[attr].apply(lambda x: "{sym}:{x}".format(sym=rep_syms[i], x=x))

        if out_hd5 is None:
            out_hd5 = copy.copy(channel)
            out_hd5["matrix"] = []
            out_hd5["barcodes"] = []

        out_df.append(df)
        out_hd5["matrix"].append(channel["matrix"])
        out_hd5["barcodes"].append(channel["barcodes"])

        print("Processed {}.".format(input_hd5_file))

    out_df = pd.concat(out_df)
    out_df.to_csv("{output_name}.attr.csv".format(output_name=output_name), na_rep="NaN")

    out_hd5["matrix"] = hstack(out_hd5["matrix"], "csc")
    out_hd5["barcodes"] = np.concatenate(out_hd5["barcodes"])
    write_10x_h5_file("{output_name}_10x.h5".format(output_name=output_name), out_hd5, genome)
    print(
        "Generated {name}_10x.h5 and {name}.attr.csv from {tot} files.".format(name=output_name, tot=len(input_names)))
