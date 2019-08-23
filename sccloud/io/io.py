#!/usr/bin/env python

import time
import numpy as np
import pandas as pd
import os.path
from scipy.io import mmread
from scipy.sparse import csr_matrix
import tables
import gzip

from typing import List, Tuple
from . import Array2D, MemData

import anndata
import logging

logger = logging.getLogger("sccloud")


def load_10x_h5_file_v2(h5_in: "tables.File", fn: str, ngene: int = None) -> "MemData":
    """Load 10x v2 format matrix from hdf5 file

    Parameters
    ----------

    h5_in : tables.File
        An instance of tables.File class that is connected to a 10x v2 formatted hdf5 file.
    fn : `str`
        File name, can be used to indicate channel-specific name prefix.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.

    Returns
    -------

    An MemData object containing genome-Array2D pair per genome.

    Examples
    --------
    >>> io.load_10x_h5_file_v2(h5_in)
    """

    data = MemData()
    for group in h5_in.list_nodes("/", "Group"):
        genome = group._v_name

        M, N = h5_in.get_node("/" + genome + "/shape").read()
        mat = csr_matrix(
            (
                h5_in.get_node("/" + genome + "/data").read(),
                h5_in.get_node("/" + genome + "/indices").read(),
                h5_in.get_node("/" + genome + "/indptr").read(),
            ),
            shape=(N, M),
        )

        barcodes = h5_in.get_node("/" + genome + "/barcodes").read().astype(str)
        ids = h5_in.get_node("/" + genome + "/genes").read().astype(str)
        names = h5_in.get_node("/" + genome + "/gene_names").read().astype(str)

        array2d = Array2D(
            {"barcodekey": barcodes}, {"featurekey": ids, "featurename": names}, mat
        )
        array2d.filter(ngene=ngene)
        array2d.separate_channels(fn)

        data.addData(genome, array2d)

    return data


def load_10x_h5_file_v3(h5_in: "tables.File", fn: str, ngene: int = None) -> "MemData":
    """Load 10x v3 format matrix from hdf5 file

    Parameters
    ----------

    h5_in : tables.File
        An instance of tables.File class that is connected to a 10x v3 formatted hdf5 file.
    fn : `str`
        File name, can be used to indicate channel-specific name prefix.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.

    Returns
    -------

    An MemData object containing genome-Array2D pair per genome.

    Examples
    --------
    >>> io.load_10x_h5_file_v3(h5_in)
    """

    M, N = h5_in.get_node("/matrix/shape").read()
    bigmat = csr_matrix(
        (
            h5_in.get_node("/matrix/data").read(),
            h5_in.get_node("/matrix/indices").read(),
            h5_in.get_node("/matrix/indptr").read(),
        ),
        shape=(N, M),
    )
    barcodes = h5_in.get_node("/matrix/barcodes").read().astype(str)
    genomes = h5_in.get_node("/matrix/features/genome").read().astype(str)
    ids = h5_in.get_node("/matrix/features/id").read().astype(str)
    names = h5_in.get_node("/matrix/features/name").read().astype(str)

    data = MemData()
    for genome in np.unique(genomes):
        idx = genomes == genome

        barcode_metadata = {"barcodekey": barcodes}
        feature_metadata = {"featurekey": ids[idx], "featurename": names[idx]}
        mat = bigmat[:, idx].copy()
        array2d = Array2D(barcode_metadata, feature_metadata, mat)
        array2d.filter(ngene)
        array2d.separate_channels(fn)

        data.addData(genome, array2d)

    return data


def load_10x_h5_file(input_h5: str, ngene: int = None) -> "MemData":
    """Load 10x format matrix (either v2 or v3) from hdf5 file

    Parameters
    ----------

    input_h5 : `str`
        The matrix in 10x v2 or v3 hdf5 format.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.

    Returns
    -------

    An MemData object containing genome-Array2D pair per genome.

    Examples
    --------
    >>> io.load_10x_h5_file('example_10x.h5')
    """

    fn = os.path.basename(input_h5)[:-3]

    data = None
    with tables.open_file(input_h5) as h5_in:
        try:
            node = h5_in.get_node("/matrix")
            data = load_10x_h5_file_v3(h5_in, fn, ngene)
        except tables.exceptions.NoSuchNodeError:
            data = load_10x_h5_file_v2(h5_in, fn, ngene)

    return data


def determine_file_name(
    path: str, names: List[str], errmsg: str, fname: str = None, exts: List[str] = None
) -> str:
    """ Try several file name options and determine which one is correct.
    """
    for name in names:
        file_name = os.path.join(path, name)
        if os.path.isfile(file_name):
            return file_name
    if fname is not None:
        for ext in exts:
            file_name = fname + ext
            if os.path.isfile(file_name):
                return file_name
    raise ValueError(errmsg)


def load_one_mtx_file(path: str, ngene: int = None, fname: str = None) -> "Array2D":
    """Load one gene-count matrix in mtx format into an Array2D object
    """
    mtx_file = determine_file_name(
        path,
        ["matrix.mtx.gz", "matrix.mtx"],
        "Expression matrix in mtx format is not found",
        fname=fname,
        exts=[".mtx"],
    )
    mat = csr_matrix(mmread(mtx_file).T)

    barcode_file = determine_file_name(
        path,
        ["cells.tsv.gz", "barcodes.tsv.gz", "barcodes.tsv"],
        "Barcode metadata information is not found",
        fname=fname,
        exts=["_barcode.tsv", ".cells.tsv"],
    )

    feature_file = determine_file_name(
        path,
        ["genes.tsv.gz", "features.tsv.gz", "genes.tsv"],
        "Feature metadata information is not found",
        fname=fname,
        exts=["_gene.tsv", ".genes.tsv"],
    )

    barcode_base = os.path.basename(barcode_file)
    feature_base = os.path.basename(feature_file)

    if barcode_base == "cells.tsv.gz" and feature_base == "genes.tsv.gz":
        format_type = "HCA DCP"
    elif barcode_base == "barcodes.tsv.gz" and feature_base == "features.tsv.gz":
        format_type = "10x v3"
    elif barcode_base == "barcodes.tsv" and feature_base == "genes.tsv":
        format_type = "10x v2"
    elif barcode_base.endswith("_barcode.tsv") and feature_base.endswith("_gene.tsv"):
        format_type = "scumi"
    elif barcode_base.endswith(".cells.tsv") and feature_base.endswith(".genes.tsv"):
        format_type = "dropEst"
    else:
        assert False

    if format_type == "HCA DCP":
        barcode_metadata = pd.read_csv(barcode_file, sep="\t", header=0)
        assert "cellkey" in barcode_metadata
        barcode_metadata.rename(columns={"cellkey": "barcodekey"}, inplace=True)

        feature_metadata = pd.read_csv(feature_file, sep="\t", header=0)
    else:
        barcode_metadata = pd.read_csv(
            barcode_file, sep="\t", header=None, names=["barcodekey"]
        )

        if format_type == "10x v3":
            feature_metadata = pd.read_csv(
                feature_file,
                sep="\t",
                header=None,
                names=["featurekey", "featurename", "featuretype"],
            )
        elif format_type == "10x v2":
            feature_metadata = pd.read_csv(
                feature_file, sep="\t", header=None, names=["featurekey", "featurename"]
            )
        elif format_type == "scumi":
            values = (
                pd.read_csv(feature_file, sep="\t", header=None)
                    .iloc[:, 0]
                    .values.astype(str)
            )
            arr = np.array(np.char.split(values, sep="_", maxsplit=1).tolist())
            feature_metadata = pd.DataFrame(
                data={"featurekey": arr[:, 0], "featurename": arr[:, 1]}
            )
        else:
            assert format_type == "dropEst"
            feature_metadata = pd.read_csv(
                feature_file, sep="\t", header=None, names=["featurekey"]
            )
            feature_metadata["featurename"] = feature_metadata["featurekey"]

    array2d = Array2D(barcode_metadata, feature_metadata, mat)
    array2d.filter(ngene=ngene)
    if format_type == "10x v3" or format_type == "10x v2":
        array2d.separate_channels("")  # fn == '' refers to 10x mtx format

    return array2d


def load_mtx_file(path: str, genome: str = None, ngene: int = None) -> "MemData":
    """Load gene-count matrix from Market Matrix files (10x v2, v3 and HCA DCP formats)

    Parameters
    ----------

    path : `str`
        Path to mtx files. The directory impiled by path should either contain matrix, feature and barcode information, or folders containg these information.
    genome : `str`, optional (default: None)
        Genome name of the matrix. If None, genome will be inferred from path.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.

    Returns
    -------

    An MemData object containing a genome-Array2D pair.

    Examples
    --------
    >>> io.load_10x_h5_file('example_10x.h5')
    """

    orig_file = None
    if not os.path.isdir(path):
        orig_file = path
        path = os.path.dirname(path)

    data = MemData()
    if (
        os.path.isfile(os.path.join(path, "matrix.mtx.gz"))
        or os.path.isfile(os.path.join(path, "matrix.mtx"))
        or (orig_file is not None and os.path.isfile(orig_file))
    ):
        if genome is None:
            genome = os.path.basename(path)
        data.addData(
            genome,
            load_one_mtx_file(
                path,
                ngene=ngene,
                fname=None if orig_file is None else os.path.splitext(orig_file)[0],
            ),
        )
    else:
        for dir_entry in os.scandir(path):
            if dir_entry.is_dir():
                data.addData(
                    dir_entry.name, load_one_mtx_file(dir_entry.path, ngene=ngene)
                )

    return data


def load_csv_file(
    input_csv: str, genome: str, sep: str = ",", ngene: int = None
) -> "MemData":
    """Load count matrix from a CSV-style file, such as CSV file or DGE style tsv file.

    Parameters
    ----------

    input_csv : `str`
        The CSV file, gzipped or not, containing the count matrix.
    genome : `str`
        The genome reference.
    sep: `str`, optional (default: ',')
        Separator between fields, either ',' or '\t'.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.

    Returns
    -------

    An MemData object containing a genome-Array2D pair.

    Examples
    --------
    >>> io.load_csv_file('example_ADT.csv', genome = 'GRCh38')
    >>> io.load_csv_file('example.umi.dge.txt.gz', genome = 'GRCh38', sep = '\t')
    """

    path = os.path.dirname(input_csv)
    base = os.path.basename(input_csv)
    is_hca_csv = base == "expression.csv"

    if sep == "\t":
        # DGE, columns are cells, which is around thousands and we can use pandas.read_csv
        df = pd.read_csv(input_csv, header=0, index_col=0, sep=sep)
        mat = csr_matrix(df.values.T)
        barcode_metadata = {"barcodekey": df.columns.values}
        feature_metadata = {
            "featurekey": df.index.values,
            "featurename": df.index.values,
        }
    else:
        # For CSV files, wide columns prevent fast pd.read_csv loading
        converter = (
            float if base.startswith("expression") else int
        )  # If expression -> float otherwise int

        barcodes = []
        names = []
        stacks = []
        with (
            gzip.open(input_csv, mode="rt")
            if input_csv.endswith(".gz")
            else open(input_csv)
        ) as fin:
            barcodes = next(fin).strip().split(sep)[1:]
            for line in fin:
                fields = line.strip().split(sep)
                names.append(fields[0])
                stacks.append([converter(x) for x in fields[1:]])

        mat = csr_matrix(np.stack(stacks, axis=1 if not is_hca_csv else 0))
        barcode_metadata = {"barcodekey": barcodes}
        feature_metadata = {"featurekey": names, "featurename": names}

        if is_hca_csv:
            barcode_file = os.path.join(path, "cells.csv")
            if os.path.exists(barcode_file):
                barcode_metadata = pd.read_csv(barcode_file, sep=",", header=0)
                assert "cellkey" in barcode_metadata
                barcode_metadata.rename(columns={"cellkey": "barcodekey"}, inplace=True)

            feature_file = os.path.join(path, "genes.csv")
            if os.path.exists(feature_file):
                feature_metadata = pd.read_csv(feature_file, sep=",", header=0)

    data = MemData()
    array2d = Array2D(barcode_metadata, feature_metadata, mat)
    array2d.filter(ngene=ngene)
    data.addData(genome, array2d)

    return data


def load_loom_file(input_loom: str, genome: str, ngene: int = None) -> "MemData":
    """Load count matrix from a LOOM file. Currently only support HCA DCP Loom spec.

    Parameters
    ----------

    input_loom : `str`
        The LOOM file, containing the count matrix.
    genome : `str`
        The genome reference.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.

    Returns
    -------

    An MemData object containing a genome-Array2D pair.

    Examples
    --------
    >>> io.load_loom_file('example.loom', genome = 'GRCh38', ngene = 200)
    """
    import loompy

    col_trans = {"CellID": "barcodekey"}
    row_trans = {"Accession": "featurekey", "Gene": "featurename"}

    data = MemData()
    with loompy.connect(input_loom) as ds:
        mat = csr_matrix(ds.sparse().T)
        barcode_metadata = {}
        for keyword, values in ds.col_attrs.items():
            keyword = col_trans.get(keyword, keyword)
            barcode_metadata[keyword] = values
        feature_metadata = {}
        for keyword, values in ds.row_attrs.items():
            keyword = row_trans.get(keyword, keyword)
            feature_metadata[keyword] = values

    array2d = Array2D(barcode_metadata, feature_metadata, mat)
    array2d.filter(ngene=ngene)
    data.addData(genome, array2d)

    return data


def load_scCloud_h5_file(
    input_h5: str, ngene: int = None, select_singlets: bool = False
) -> "MemData":
    """Load matrices from sccloud-format hdf5 file

    Parameters
    ----------

    input_h5 : `str`
        sccloud-format hdf5 file.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.
    select_singlets: `bool`, optional (default: False)
        If only load singlets.

    Returns
    -------

    An MemData object containing genome-Array2D pair per genome.

    Examples
    --------
    >>> io.load_scCloud_h5_file('example.sccloud.h5')
    """

    cite_seq_name = None
    selected_barcodes = None

    data = MemData()
    with tables.open_file(input_h5) as h5_in:
        for group in h5_in.list_nodes("/", "Group"):
            genome = group._v_name

            M, N = h5_in.get_node("/" + genome + "/shape").read()
            mat = csr_matrix(
                (
                    h5_in.get_node("/" + genome + "/data").read(),
                    h5_in.get_node("/" + genome + "/indices").read(),
                    h5_in.get_node("/" + genome + "/indptr").read(),
                ),
                shape=(N, M),
            )

            barcode_metadata = {}
            for node in h5_in.walk_nodes("/" + genome + "/_barcodes", "Array"):
                values = node.read()
                if values.dtype.kind == "S":
                    values = values.astype(str)
                barcode_metadata[node.name] = values

            feature_metadata = {}
            for node in h5_in.walk_nodes("/" + genome + "/_features", "Array"):
                values = node.read()
                if values.dtype.kind == "S":
                    values = values.astype(str)
                feature_metadata[node.name] = values

            array2d = Array2D(barcode_metadata, feature_metadata, mat)
            if genome.startswith("CITE_Seq"):
                cite_seq_name = genome
            else:
                array2d.filter(ngene, select_singlets)
                selected_barcodes = array2d.get_metadata("barcodekey")

            data.addData(genome, array2d)

    if (cite_seq_name is not None) and (selected_barcodes is not None):
        array2d = data.getData(cite_seq_name)
        selected = array2d.get_metadata("barcodekey").isin(selected_barcodes)
        array2d.trim(selected)

    return data


def infer_file_format(input_file: str) -> Tuple[str, str, str]:
    """ Infer file format from input_file name

    This function infer file format by inspecting the file name.

    Parameters
    ----------

    input_file : `str`
        Input file name.

    Returns
    -------
    `str`
        File format, choosing from 'sccloud', '10x', 'h5ad', 'mtx', 'dge', and 'csv'.
    `str`
        The path covering all input files. Most time this is the same as input_file. But for HCA mtx and csv, this should be parent directory.
    `str`
        Type of the path, either 'file' or 'directory'.
    """

    file_format = None
    copy_path = input_file
    copy_type = "file"

    if input_file.endswith(".h5"):
        file_format = "10x"
    elif input_file.endswith(".h5sc"):
        file_format = "sccloud"
    elif input_file.endswith(".h5ad"):
        file_format = "h5ad"
    elif input_file.endswith(".loom"):
        file_format = "loom"
    elif (
        input_file.endswith(".mtx")
        or input_file.endswith(".mtx.gz")
        or os.path.splitext(input_file)[1] == ""
    ):
        file_format = "mtx"
        if os.path.splitext(input_file)[1] != "":
            copy_path = os.path.dirname(input_file)
        copy_type = "directory"
    elif input_file.endswith("dge.txt.gz"):
        file_format = "dge"
    elif input_file.endswith(".csv") or input_file.endswith(".csv.gz"):
        file_format = "csv"
        if os.path.basename(input_file) == "expression.csv":
            copy_type = os.path.dirname(input_file)
            copy_type = "directory"
    else:
        raise ValueError("Unrecognized file type for file {}!".format(input_file))

    return file_format, copy_path, copy_type


def read_input(
    input_file: str,
    genome: str = None,
    return_type: str = "AnnData",
    concat_matrices: bool = False,
    h5ad_mode: str = "a",
    ngene: int = None,
    select_singlets: bool = False,
) -> "MemData or AnnData or List[AnnData]":
    """Load data into memory.

    This function is used to load input data into memory. Inputs can be in 10x genomics v2 & v3 formats (hdf5 or mtx), HCA DCP mtx and csv formats, Drop-seq dge format, and CSV format.

    Parameters
    ----------

    input_file : `str`
        Input file name.
    genome : `str`, optional (default: None)
        A string contains comma-separated genome names. sccloud will read all matrices matching the genome names. If genomes is None, all matrices will be considered.
    return_type : `str`
        Return object type, can be either 'MemData' or 'AnnData'.
    concat_matrices : `boolean`, optional (default: False)
        If input file contains multiple matrices, if concatenate them into one AnnData object or return a list of AnnData objects.
    h5ad_mode : `str`, optional (default: `a`)
        If input is in h5ad format, the backed mode for loading the data. mode could be 'a', 'r', 'r+'. 'a' refers to load all into memory.
    ngene : `int`, optional (default: None)
        Minimum number of genes to keep a barcode. Default is to keep all barcodes.
    select_singlets : `bool`, optional (default: False)
        If only keep DemuxEM-predicted singlets when loading data.

    Returns
    -------
    `MemData` object or `anndata` object or a list of `anndata` objects
        An `MemData` object or `anndata` object or a list of `anndata` objects containing the count matrices.

    Examples
    --------
    >>> adata = io.read_input('example_10x.h5', genomes = 'mm10')
    >>> adata = io.read_input('example.h5ad', mode = 'r+')
    >>> adata = io.read_input('example_ADT.csv')
    """

    start = time.time()

    input_file = os.path.expanduser(os.path.expandvars(input_file))
    file_format, _, _ = infer_file_format(input_file)

    if file_format == "sccloud":
        data = load_scCloud_h5_file(
            input_file, ngene=ngene, select_singlets=select_singlets
        )
    elif file_format == "10x":
        data = load_10x_h5_file(input_file, ngene=ngene)
    elif file_format == "h5ad":
        data = anndata.read_h5ad(
            input_file, backed=(None if h5ad_mode == "a" else h5ad_mode)
        )
    elif file_format == "mtx":
        data = load_mtx_file(input_file, genome, ngene=ngene)
    elif file_format == "loom":
        assert genome is not None
        data = load_loom_file(input_file, genome, ngene=ngene)
    else:
        assert (file_format == "dge" or file_format == "csv") and (genome is not None)
        data = load_csv_file(
            input_file, genome, sep=("\t" if file_format == "dge" else ","), ngene=ngene
        )

    if file_format != "h5ad":
        data.restrain_keywords(genome)
        if return_type == "AnnData":
            data = data.convert_to_anndata(concat_matrices=concat_matrices)
    else:
        assert return_type == "AnnData"

    end = time.time()
    logger.info("Read input is finished. Time spent = {:.2f}s.".format(end - start))

    return data


def write_output(data: "MemData or AnnData", output_name: str) -> None:
    """ Write data back to disk.

    This function is used to write data back to disk.

    Parameters
    ----------

    data : `MemData` or `AnnData`
        data to write back, can be either an MemData or AnnData object.
    output_name : `str`
        output file name. MemData ends with suffix '.h5sc' and AnnData ends with suffix '.h5ad'. If output_name has no suffix, an appropriate suffix will be appended.

    Returns
    -------
    None

    Examples
    --------
    >>> io.write_output(adata, 'test')
    """

    start = time.time()

    if isinstance(data, MemData):
        if not output_name.endswith(".h5sc"):
            output_name += ".h5sc"
        data.write_h5_file(output_name)
    else:
        output_file_format = 'h5ad'
        output_name_lc = output_name.lower()
        file_formats = ['h5ad', 'loom']
        for file_format in file_formats:
            if output_name_lc.endswith('.' + file_format):
                output_file_format = file_format
                break

        if not output_name_lc.endswith("." + output_file_format):
            output_name += "." + output_file_format

        # Eliminate objects starting with fmat_ from uns
        keys = list(data.uns)
        for keyword in keys:
            if keyword.startswith("fmat_"):
                data.uns.pop(keyword)

        if output_file_format == 'h5ad':
            import pathlib

            output_name = pathlib.Path(output_name)
            if data.isbacked and output_name == data.filename:
                import h5py

                h5_file = data.file._file
                # Fix old h5ad files in which obsm/varm were stored as compound datasets
                if "obsm" in h5_file.keys() and isinstance(h5_file["obsm"], h5py.Dataset):
                    del h5_file["obsm"]
                if "varm" in h5_file.keys() and isinstance(h5_file["varm"], h5py.Dataset):
                    del h5_file["varm"]
            data.write(output_name, compression="gzip")
        elif output_file_format == 'loom':
            data.write_loom(output_name, write_obsm_varm=True)
        else:
            raise ValueError('Unknown file format')

    end = time.time()
    logger.info("Write output is finished. Time spent = {:.2f}s.".format(end - start))
