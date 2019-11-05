#!/usr/bin/env python

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, hstack, vstack
import tables

from typing import List
from pegasus.utils import decorators as pg_deco

import anndata


class Array2D:
    def __init__(
        self,
        barcode_metadata: "dict or pd.DataFrame" = {},
        feature_metadata: "dict or pd.DataFrame" = {},
        matrix: "csr_matrix" = csr_matrix((0, 0)),
        metadata: dict = {},
    ):
        self.barcode_metadata = (
            barcode_metadata
            if isinstance(barcode_metadata, pd.DataFrame)
            else pd.DataFrame(barcode_metadata)
        )
        self.feature_metadata = (
            feature_metadata
            if isinstance(feature_metadata, pd.DataFrame)
            else pd.DataFrame(feature_metadata)
        )
        self.matrix = matrix  # scipy csr matrix
        self.metadata = metadata  # other metadata, a dictionary

        if "barcodekey" in self.barcode_metadata:
            self.barcode_metadata.set_index("barcodekey", inplace=True)
        else:
            assert self.barcode_metadata.index.name == "barcodekey"
        if "featurekey" in self.feature_metadata:
            self.feature_metadata.set_index("featurekey", inplace=True)
        else:
            assert self.feature_metadata.index.name == "featurekey"

        if self.barcode_metadata.shape[0] != self.matrix.shape[0]:
            raise ValueError(
                "Wrong number of cells : matrix has {} cells, barcodes file has {}".format(
                    self.matrix.shape[0], self.barcode_metadata.shape[0]
                )
            )
        if self.feature_metadata.shape[0] != (
            self.matrix.shape[1] if len(self.matrix.shape) == 2 else 1
        ):
            raise ValueError(
                "Wrong number of features : matrix has {} features, features file has {}".format(
                    self.matrix.shape[1], self.feature_metadata.shape[0]
                )
            )

    def list_metadata_keys(self, type: str = "barcode") -> List[str]:
        """ Return available keys in metadata, type = barcode or feature or None
        """
        if type == "barcode":
            return [
                self.barcode_metadata.index.name
            ] + self.barcode_metadata.columns.tolist()
        elif type == "feature":
            return [
                self.feature_metadata.index.name
            ] + self.feature_metadata.columns.tolist()
        else:
            list(self.metadata)

    def get_metadata(self, name: str, type: str = "barcode") -> List[object]:
        """ Return values of column name in metadata type, type can be barcode, feature, or None
        """
        if type == "barcode":
            return (
                self.barcode_metadata.index
                if name == "barcodekey"
                else self.barcode_metadata[name].values
            )
        elif type == "feature":
            return (
                self.feature_metadata.index
                if name == "featurekey"
                else self.feature_metadata[name].values
            )
        else:
            return self.metadata[name]

    def trim(self, selected: List[bool]) -> None:
        """ Only keep barcodes in selected
        """
        self.matrix = self.matrix[selected, :]
        self.barcode_metadata = self.barcode_metadata[selected]

    def filter(self, ngene: int = None, select_singlets: bool = False) -> None:
        """ Filter out low quality barcodes, only keep barcodes satisfying ngene >= ngene and selecting singlets if select_singlets is True
        """
        if (ngene is None) and (not select_singlets):
            return None

        selected = np.ones(self.matrix.shape[0], dtype=bool)
        if ngene is not None:
            selected = selected & (self.matrix.getnnz(axis=1) >= ngene)
        if select_singlets:
            assert "demux_type" in self.barcode_metadata
            selected = (
                selected & (self.barcode_metadata["demux_type"] == "singlet").values
            )
            self.barcode_metadata.drop(columns="demux_type", inplace=True)

        self.trim(selected)

    def separate_channels(self, fn: str) -> None:
        """ Separate channel information from barcodekeys, only used for 10x v2, v3 h5 and mtx.
        """

        barcodes = self.barcode_metadata.index.values.astype(str)
        parts = np.array(np.char.rsplit(barcodes, sep="-", maxsplit=1).tolist())

        channels = None
        assert fn == "" or parts.shape[1] == 2
        if fn != "" and np.char.not_equal(parts[:, 1], "1").sum() > 0:
            # if we have multiple channels
            channels = [fn + "-" + x for x in parts[:, 1]]
            barcodes = np.array([x + "-" + y for x, y in zip(channels, parts[:, 0])])
        else:
            barcodes = parts[:, 0]

        self.barcode_metadata.index = pd.Index(barcodes, name="barcodekey")

        if channels is not None:
            self.barcode_metadata["Channel"] = channels

    def update_barcode_metadata_info(
        self, sample_name: str, row: "pd.Series", attributes: List[str]
    ) -> None:
        """ Update barcodekey, update channel and add attributes
        """
        nsample = self.barcode_metadata.shape[0]
        barcodes = [sample_name + "-" + x for x in self.barcode_metadata.index]
        self.barcode_metadata.index = pd.Index(barcodes, name="barcodekey")
        if "Channel" in self.barcode_metadata:
            self.barcode_metadata["Channel"] = [
                sample_name + "-" + x for x in self.barcode_metadata["Channel"]
            ]
        else:
            self.barcode_metadata["Channel"] = np.repeat(sample_name, nsample)
        if attributes is not None:
            for attr in attributes:
                self.barcode_metadata[attr] = np.repeat(row[attr], nsample)

    def write_to_hdf5(self, keyword: str, hd5_out: "File", is_h5sc: bool) -> None:
        """ Write Array2D content into hdf5 file
        """
        out_group = hd5_out.create_group("/", keyword)
        # write matrix
        hd5_out.create_carray(out_group, "data", obj=self.matrix.data)
        hd5_out.create_carray(out_group, "indices", obj=self.matrix.indices)
        hd5_out.create_carray(out_group, "indptr", obj=self.matrix.indptr)
        M, N = self.matrix.shape
        hd5_out.create_carray(
            out_group, "shape", obj=(N, M)
        )  # store as feature by barcode instead
        # write barcode_metadata
        outgb = hd5_out.create_group("/" + keyword, "_barcodes") if is_h5sc else out_group
        if not is_h5sc:
            self.barcode_metadata.index.name = "barcodes"
            self.feature_metadata.index.name = "genes"
            self.feature_metadata.rename(columns = {"featurename" : "gene_names"}, inplace = True)
        hd5_out.create_carray(
            outgb,
            self.barcode_metadata.index.name,
            obj=self.barcode_metadata.index.values.astype("S"),
        )  # encode into binary strings
        for col in self.barcode_metadata:
            kind = self.barcode_metadata[col].dtype.kind
            if kind == "U" or kind == "O":
                hd5_out.create_carray(
                    outgb, col, obj=self.barcode_metadata[col].values.astype("S")
                )  # encode into binary strings
            else:
                hd5_out.create_carray(outgb, col, obj=self.barcode_metadata[col].values)
        # write feature_metadata
        outgb = hd5_out.create_group("/" + keyword, "_features") if is_h5sc else out_group
        hd5_out.create_carray(
            outgb,
            self.feature_metadata.index.name,
            obj=self.feature_metadata.index.values.astype("S"),
        )  # encode into binary strings
        for col in self.feature_metadata:
            kind = self.feature_metadata[col].dtype.kind
            if kind == "U" or kind == "O":
                hd5_out.create_carray(
                    outgb, col, obj=self.feature_metadata[col].values.astype("S")
                )  # encode into binary strings
            else:
                hd5_out.create_carray(outgb, col, obj=self.feature_metadata[col].values)


def polish_featurename(
    feature_names: List[str], feature_keys: List[str], genomes: List[str]
) -> List[str]:
    """    Remove prefixing genome strings and deduplicate feature names
    """
    import re
    from collections import Counter

    prefix = re.compile("^(" + "|".join(genomes) + ")_+")
    if prefix.match(feature_names[0]):
        feature_names = [prefix.sub("", x) for x in feature_names]
        feature_keys = [prefix.sub("", x) for x in feature_keys]

    dup_ids = Counter()
    for i in range(len(feature_names)):
        idn = dup_ids[feature_names[i]]
        dup_ids[feature_names[i]] += 1
        if idn > 0:
            feature_names[i] = feature_names[i] + ".#~{}".format(idn + 1) # duplicate ID starts from 2. .#~ makes it unique.

    return feature_names, feature_keys


def get_fillna_dict(df: "pd.DataFrame") -> dict:
    """ Generate a fillna dict for columns in a df
    """
    fillna_dict = {}
    for column in df:
        if df[column].dtype.kind in {"O", "S"}:
            fillna_dict[column] = ""
        else:
            fillna_dict[column] = 0
    return fillna_dict


class MemData:
    def __init__(self):
        self.data = {}  # data is a dictionary mapping keyword to Array2D

    def listKeys(self) -> List[str]:
        return list(self.data)

    def getData(self, keyword: str) -> "Array2D":
        assert keyword in self.data
        return self.data[keyword]

    def addData(self, keyword: str, one_array: "Array2D") -> None:
        self.data[keyword] = one_array

    def update_barcode_metadata_info(
        self, sample_name: str, row: "pd.Series", attributes: List[str]
    ) -> None:
        """ Update barcodekey, update channel and add attributes for each array2d array
        """
        for array2d in self.data.values():
            array2d.update_barcode_metadata_info(sample_name, row, attributes)

    def addAggrData(self, data: "MemData") -> None:
        """ Add Aggr Data
        """
        for keyword, array2d in data.data.items():
            if keyword in self.data:
                self.data[keyword].append(array2d)
            else:
                self.data[keyword] = [array2d]

    @pg_deco.TimeLogger()
    def aggregate(self) -> None:
        """ Merge aggregated count matrices
        """
        import gc
        import warnings

        for keyword in self.data:
            array2d_list = self.data[keyword]
            self.data[keyword] = None

            if len(array2d_list) == 1:
                self.data[keyword] = array2d_list[0]
            else:
                barcode_metadata_dfs = [
                    array2d.barcode_metadata for array2d in array2d_list
                ]
                barcode_metadata = pd.concat(barcode_metadata_dfs, axis=0, sort=False)
                fillna_dict = get_fillna_dict(barcode_metadata)
                barcode_metadata.fillna(value=fillna_dict, inplace=True)

                feature_metadata = array2d_list[0].feature_metadata
                for other in array2d_list[1:]:
                    keys = ["featurekey"] + feature_metadata.columns.intersection(
                        other.feature_metadata.columns
                    ).values.tolist()
                    feature_metadata = feature_metadata.merge(
                        other.feature_metadata, on=keys, how="outer", sort=False
                    )  # If sort is True, feature keys will be changed even if all channels share the same feature keys.
                fillna_dict = get_fillna_dict(feature_metadata)
                feature_metadata.fillna(value=fillna_dict, inplace=True)

                matrix_list = []
                f2idx = pd.Series(
                    data=range(feature_metadata.shape[0]), index=feature_metadata.index
                )
                for array2d in array2d_list:
                    if (
                        feature_metadata.shape[0] > array2d.feature_metadata.shape[0]
                        or (
                            feature_metadata.index != array2d.feature_metadata.index
                        ).sum()
                        > 0
                    ):
                        mat = csr_matrix(
                            (array2d.matrix.shape[0], f2idx.size),
                            dtype=array2d.matrix.dtype,
                        )

                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            mat[
                                :, f2idx[array2d.feature_metadata.index].values
                            ] = array2d.matrix

                        array2d.matrix = mat
                        gc.collect()

                    matrix_list.append(array2d.matrix)

                newmat = vstack(matrix_list)
                matrix_list = array2d_list = None
                gc.collect()
                self.data[keyword] = Array2D(barcode_metadata, feature_metadata, newmat)

    def restrain_keywords(self, keywords: str) -> None:
        if keywords is None:
            return None

        keywords = set(keywords.split(","))
        available = set(self.data)

        invalid_set = keywords - available
        if len(invalid_set) > 0:
            raise ValueError(
                "Keywords {} do not exist.".format(",".join(list(invalid_set)))
            )

        remove_set = available - keywords
        for keyword in remove_set:
            self.data.pop(keyword)

    def convert_to_anndata(self, concat_matrices: bool =False, channel_attr: str = None, black_list: List[str] = []) -> "AnnData or List[AnnData]":
        """Convert an MemData object into SCANPY's AnnData object

        Parameters
        ----------

        concat_matrices : `bool`, optional (default: False)
            If concatenate multiple matrices. If so, return only one AnnData object, otherwise, might return a list of AnnData objects.

        channel_attr : `str`, optional (default: None)
            Use channel_attr to represent different samples. This will set a 'Channel' column field with channel_attr.

        black_list : `List[str]`, optional (default: [])
            Attributes in black list will be poped out.

        Returns
        -------

        `anndata` object or a dictionary of `anndata` objects
            An `anndata` object or a list of `anndata` objects containing the count matrices.

        Warning
        -------
        This procedure will convert all int matrix into float32 matrix in place!

        Examples
        --------
        >>> adata = MemData.convert_to_anndata()
        """

        Xs = []
        feature_dfs = []
        results = []

        genomes = self.listKeys()
        for genome in genomes:
            array2d = self.data[genome]

            if (
                array2d.matrix.dtype == np.int32
            ):  # caution, change matrix from int to float32 permanently!
                array2d.matrix.dtype = np.float32
                orig_data = array2d.matrix.data.view(np.int32)
                array2d.matrix.data[:] = orig_data

            obs_dict = array2d.barcode_metadata.to_dict(orient="list")
            obs_dict["obs_names"] = array2d.barcode_metadata.index.values
            if (channel_attr is not None) and (channel_attr in obs_dict):
                obs_dict["Channel"] = obs_dict[channel_attr]
            if "Channel" not in obs_dict:
                obs_dict["Channel"] = [""] * array2d.barcode_metadata.shape[0]
            for attr in black_list:
                if attr in obs_dict:
                    obs_dict.pop(attr)

            if not concat_matrices or len(genomes) == 1:
                var_dict = array2d.feature_metadata.to_dict(orient="list")
                feature_names, feature_keys = polish_featurename(
                    var_dict.pop("featurename"),
                    array2d.feature_metadata.index.values,
                    [genome],
                )
                var_dict["var_names"] = feature_names
                var_dict["gene_ids"] = feature_keys

                adata = anndata.AnnData(X=array2d.matrix, obs=obs_dict, var=var_dict)
                adata.uns["genome"] = genome
                results.append(adata)
            else:
                Xs.append(array2d.matrix)
                feature_dfs.append(array2d.feature_metadata)

        if len(results) == 0:
            fdf = pd.concat(feature_dfs, axis=0)
            fdf.fillna(value="N/A", inplace=True)
            var_dict = fdf.to_dict(orient="list")
            feature_names, feature_keys = polish_featurename(
                var_dict.pop("featurename"), fdf.index.values, genomes
            )
            var_dict["var_names"] = feature_names
            var_dict["gene_ids"] = feature_keys

            results = anndata.AnnData(
                X=hstack(Xs, format="csr"), obs=obs_dict, var=var_dict
            )
            results.uns["genome"] = ",".join(genomes)
        elif len(results) == 1:
            results = results[0]

        return results

    def write_h5_file(self, output_h5: str) -> None:
        """Write count matricies into a pegasus-format HDF5 file --- each matrix is stored in a group under root. For each matrix, the barcode metadata are stored under _barcodes subgroup and the feature metadata are stored under _features subgroup

        Parameters
        ----------

        output_h5 : `str`
            The output file name. If ends with .h5sc, output h5sc format; if ends with .h5, output 10x V2 h5 format.

        Returns
        -------

        None

        Examples
        --------
        >>> MemData.write_h5_file('example_10x.h5')
        """

        is_h5sc = output_h5.endswith(".h5sc")
        with tables.open_file(
            output_h5, mode="w", title=output_h5, filters=tables.Filters(complevel=1)
        ) as hd5_out:
            for keyword, array2d in self.data.items():
                array2d.write_to_hdf5(keyword, hd5_out, is_h5sc)

