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


from pegasusio import infer_file_type, read_input, write_output
from .tools import (
    qc_metrics,
    get_filter_stats,
    filter_data,
    identify_robust_genes,
    log_norm,
    select_features,
    pca,
    highly_variable_features,
    set_group_attribute,
    correct_batch,
    run_harmony,
    run_scanorama,
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
from .misc import search_genes, search_de_genes
from .plotting import (
    scatter,
    compo_plot,
    scatter_groups,
    violin,
    heatmap,
    dotplot,
    dendrogram,
    hvfplot,
    qcviolin,
    volcano,
)

from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version('pegasuspy')
    del version
except PackageNotFoundError:
    pass
