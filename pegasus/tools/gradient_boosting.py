import time
import numpy as np
import pandas as pd
from collections import defaultdict
from joblib import effective_n_jobs

from typing import List, Dict
from anndata import AnnData

from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans

from pegasusio import read_input

import logging
logger = logging.getLogger(__name__)

from pegasusio import timer


@timer(logger=logger)
def find_markers(
    data: AnnData,
    label_attr: str,
    de_key: str = "de_res",
    n_jobs: int = -1,
    min_gain: float = 1.0,
    random_state: int = 0,
    remove_ribo: bool = False,
) -> Dict[str, Dict[str, List[str]]]:
    """Find markers using gradient boosting method.

    Parameters
    ----------
    data: ``anndata.AnnData``
        Annotated data matrix with rows for cells and columns for genes.

    label_attr: ``str``
        Cluster labels used for finding markers. Must exist in ``data.obs``.

    de_key: ``str``, optional, default: ``"de_res"``
        Keyword of DE analysis result stored in ``data.varm``.

    n_jobs: ``int``, optional, default: ``-1``
        Number of threads to used. If ``-1``, use all available threads.

    min_gain: ``float``, optional, default: ``1.0``
        Only report genes with a feature importance score (in gain) of at least ``min_gain``.

    random_state: ``int``, optional, default: ``0``
        Random seed set for reproducing results.

    remove_ribo: ``bool``, optional, default: ``False``
        If ``True``, remove ribosomal genes with either RPL or RPS as prefixes.

    Returns
    -------
    markers: ``Dict[str, Dict[str, List[str]]]``
        A Python dictionary containing marker information in structure ``dict[cluster_id]['up' or 'down'][dataframe]``.

    Examples
    --------
    >>> marker_dict = pg.find_markers(adata, label_attr = 'leiden_labels')
    """

    n_jobs = effective_n_jobs(n_jobs)

    if remove_ribo:
        data = data[
            :,
            np.vectorize(lambda x: not x.startswith("RPL") and not x.startswith("RPS"))(
                data.var_names
            ),
        ]

    X_train, X_test, y_train, y_test = train_test_split(
        data.X,
        data.obs[label_attr],
        test_size=0.1,
        random_state=random_state,
        stratify=data.obs[label_attr],
    )

    # start = time.time()
    # xgb = XGBClassifier(n_jobs = n_jobs, n_gpus = 0)
    # xgb.fit(X_train, y_train, eval_set = [(X_train, y_train), (X_test, y_test)], eval_metric = 'merror')
    # # print(xgb.evals_result())
    # end = time.time()
    # print("XGBoost used {:.2f}s to train.".format(end - start))

    # from xgboost import XGBClassifier
    try:
        from lightgbm import LGBMClassifier
    except ImportError:
        print("Need lightgbm! Try 'pip install lightgbm'.")
    start_lgb = time.time()
    lgb = LGBMClassifier(n_jobs=n_jobs, metric="multi_error", importance_type="gain")
    lgb.fit(
        X_train,
        y_train,
        eval_set=[(X_train, y_train), (X_test, y_test)],
        early_stopping_rounds=1,
    )
    end_lgb = time.time()
    logger.info("LightGBM used {:.2f}s to train.".format(end_lgb - start_lgb))

    ntot = (lgb.feature_importances_ >= min_gain).sum()
    ords = np.argsort(lgb.feature_importances_)[::-1][:ntot]

    log_exprs = [
        x for x in data.varm[de_key].dtype.names if x.startswith("mean_logExpr:")
    ]
    labels = [x.rpartition(":")[2] for x in log_exprs]

    titles = [("down", "down_gain"), ("weak", "weak_gain"), ("strong", "strong_gain")]
    markers = defaultdict(lambda: defaultdict(list))

    kmeans = KMeans(n_clusters=3, random_state=random_state)
    for gene_id in ords:
        gene_symbol = data.var_names[gene_id]
        mydat = [[x] for x in data.varm[de_key][log_exprs][gene_id]]
        kmeans.fit(mydat)
        kmeans_label_mode = pd.Series(kmeans.labels_).mode()[0]
        for i, kmeans_label in enumerate(np.argsort(kmeans.cluster_centers_[:, 0])):
            if kmeans_label != kmeans_label_mode:
                for pos in (kmeans.labels_ == kmeans_label).nonzero()[0]:
                    clust_label = labels[pos]
                    markers[clust_label][titles[i][0]].append(gene_symbol)
                    markers[clust_label][titles[i][1]].append(
                        "{:.2f}".format(lgb.feature_importances_[gene_id])
                    )

    return markers


def run_find_markers(
    input_h5ad_file: str,
    output_file: str,
    label_attr: str,
    de_key: str = "de_res",
    n_jobs: int = -1,
    min_gain: float = 1.0,
    random_state: int = 0,
    remove_ribo: bool = False,
) -> None:
    """
    For command line use.
    """
    import xlsxwriter
    from natsort import natsorted

    data = read_input(input_h5ad_file)
    markers = find_markers(
        data,
        label_attr,
        de_key=de_key,
        n_jobs=n_jobs,
        min_gain=min_gain,
        random_state=random_state,
        remove_ribo=remove_ribo,
    )

    keywords = [("strong", "strong_gain"), ("weak", "weak_gain"), ("down", "down_gain")]

    writer = pd.ExcelWriter(output_file, engine="xlsxwriter")

    for clust_id in natsorted(markers.keys()):
        clust_markers = markers[clust_id]

        sizes = []
        for keyword in keywords:
            sizes.append(len(clust_markers[keyword[0]]))

        arr = np.zeros((max(sizes), 8), dtype=object)
        arr[:] = ""

        for i in range(3):
            arr[0 : sizes[i], i * 3] = clust_markers[keywords[i][0]]
            arr[0 : sizes[i], i * 3 + 1] = clust_markers[keywords[i][1]]

        df = pd.DataFrame(
            data=arr,
            columns=[
                "strongly up-regulated",
                "gain",
                "",
                "weakly up-regulated",
                "gain",
                "",
                "down-regulated",
                "gain",
            ],
        )
        df.to_excel(writer, sheet_name=clust_id, index=False)

    writer.save()
