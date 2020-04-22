import logging
from typing import Dict

import pandas as pd
import time

logger = logging.getLogger("pegasus")


def enrichment_analysis(markers: Dict[str, Dict[str, pd.DataFrame]], max_genes: int = 100, organism: str = 'hsapiens',
                        enrichment_threshold: float = 0.05) -> pd.DataFrame:
    """Perform enrichment analysis using gprofiler (https://biit.cs.ut.ee/gprofiler/gost).

    Parameters
    ----------
    markers: ``Dict[str, Dict[str, pd.DataFrame]``
        Output from markers.

    max_genes: ``int``, optional, default: 100
        Maximum number of genes to use in enrichment query
    organism: ``str``, optional, default: ``hsapiens``
        Organism. See https://biit.cs.ut.ee/gprofiler/page/organism-list for full list.
    enrichment_threshold: ``float``, optional, default: ``0.05``
        Include enrichment results with corrected p-value less than this threshold

    Returns
    -------
    ``pd.DataFrame``

    """
    start = time.perf_counter()
    from gprofiler import GProfiler
    gp = GProfiler(return_dataframe=True)
    query = {}
    for cluster in markers.keys():
        up_list = markers[cluster]['up'].index.values.tolist()
        if len(up_list) > 0:
            query[cluster + '-up'] = up_list[0:max_genes]
        down_list = markers[cluster]['down'].index.values.tolist()
        if len(down_list) > 0:
            query[cluster + '-down'] = down_list[0:max_genes]
    result = gp.profile(organism=organism, query=query, user_threshold=enrichment_threshold)
    end = time.perf_counter()
    logger.info(
        "Enrichment analysis is finished. Time spent = {:.2f}s.".format(
            end - start
        )
    )
    return result
