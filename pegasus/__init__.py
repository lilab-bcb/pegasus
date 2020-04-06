try:
    get_ipython
except NameError:
    import matplotlib

    matplotlib.use("Agg")

import sys
import logging
import warnings

logger = logging.getLogger("pegasus")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

warnings.filterwarnings("ignore", category=UserWarning,  module='lightgbm')
warnings.filterwarnings("ignore", category=FutureWarning, module='anndata')


from .io import infer_file_format, read_input, write_output
from .tools import (
    aggregate_matrices,
    qc_metrics,
    get_filter_stats,
    filter_data,
    log_norm,
    select_features,
    pca,
    highly_variable_features,
    set_group_attribute,
    correct_batch,
    run_harmony,
    neighbors,
    calc_kBET,
    calc_kSIM,
    diffmap,
    reduce_diffmap_to_3d,
    calc_pseudotime,
    cluster,
    louvain,
    leiden,
    spectral_louvain,
    spectral_leiden,
    tsne,
    fitsne,
    umap,
    fle,
    net_tsne,
    net_fitsne,
    net_umap,
    net_fle,
    de_analysis,
    markers,
    write_results_to_excel,
    find_markers,
    infer_path,
    calc_signature_score,
)
from .annotate_cluster import infer_cell_types, annotate, infer_cluster_names
from .demuxEM import estimate_background_probs, demultiplex
from .misc import search_genes, search_de_genes

from scplot import (
    violin,
    heatmap,
    scatter,
    line,
    dotplot,
    scatter_matrix,
    embedding,
    composition_plot,
    variable_feature_plot,
    volcano,
)

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
