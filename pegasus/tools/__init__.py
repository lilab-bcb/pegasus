import psutil
_cpu_count = psutil.cpu_count(logical=False)
if _cpu_count is None:
    _cpu_count = psutil.cpu_count(logical=True)


from .utils import eff_n_jobs, update_rep, X_from_rep, W_from_rep, slicing, calc_mean, calc_mean_and_var, calc_expm1, calc_stat_per_batch, normalize_by_count, calc_sig_background, simulate_doublets, check_batch_key

from .preprocessing import (
    qc_metrics,
    get_filter_stats,
    filter_data,
    identify_robust_genes,
    _run_filter_data,
    log_norm,
    select_features,
    pca,
    pc_transform,
    tsvd,
    tsvd_transform,
    regress_out,
)
from .hvf_selection import estimate_feature_statistics, highly_variable_features
from .batch_correction import run_harmony, run_scanorama
from .nearest_neighbors import (
    calculate_nearest_neighbors,
    get_neighbors,
    neighbors,
    calculate_affinity_matrix,
    calc_kBET,
    calc_kSIM,
)
from .graph_operations import construct_graph
from .diffusion_map import diffmap
from .pseudotime import calc_pseudotime, infer_path
from .clustering import jump_method, louvain, leiden, spectral_louvain, spectral_leiden, cluster, split_one_cluster
from .net_regressor import net_train_and_predict
from .visualization import (
    tsne,
    umap,
    fle,
    net_umap,
    net_fle,
)
from .diff_expr import de_analysis, markers, write_results_to_excel, run_de_analysis
from .gradient_boosting import find_markers, run_find_markers
from .subcluster_utils import clone_subset
from .signature_score import calc_signature_score, calculate_z_score
from .doublet_detection import infer_doublets, mark_doublets
from .nmf import nmf, integrative_nmf
