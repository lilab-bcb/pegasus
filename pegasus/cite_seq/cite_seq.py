import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from pegasus.io import read_input, write_output


def merge_rna_and_adt_data(input_raw_h5, input_csv, antibody_control_csv, output_name):
    data = read_input(input_raw_h5, return_type="MemData")
    print("Loaded the RNA matrix.")

    keyword = "CITE_Seq_" + data.listKeys()[0]
    data_citeseq = read_input(input_csv, return_type="MemData", genome=keyword)
    print("Loaded the ADT matrix.")

    array2d = data_citeseq.getData(keyword)
    if antibody_control_csv is None:
        array2d.matrix = array2d.matrix.log1p()
    else:
        size = array2d.feature_metadata.shape[0]
        idx = np.zeros(size, dtype=bool)
        antibody_to_pos = pd.Series(
            data=range(size), index=array2d.feature_metadata.index
        )

        adt_matrix = array2d.matrix.toarray().astype(float)

        series = pd.read_csv(antibody_control_csv, header=0, index_col=0, squeeze=True)
        for antibody, control in series.iteritems():
            pos_a = antibody_to_pos[antibody]
            pos_c = antibody_to_pos[control]
            idx[pos_a] = True
            # convert to log expression
            adt_matrix[:, pos_a] = np.maximum(
                np.log(adt_matrix[:, pos_a] + 1.0) - np.log(adt_matrix[:, pos_c] + 1.0),
                0.0,
            )

        array2d.feature_metadata = array2d.feature_metadata[idx]
        array2d.matrix = csr_matrix(adt_matrix[:, idx])

    data.addData(keyword, array2d)
    write_output(data, output_name)

    print("Merged output is written.")


def capping(adt_matrix, percentile):
    for i in range(adt_matrix.shape[1]):
        cap = np.percentile(adt_matrix[:, i], percentile)
        adt_matrix[adt_matrix[:, i] > cap, i] = cap
