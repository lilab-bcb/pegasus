import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix, hstack

from pegasusio import UnimodalData, MultimodalData
from pegasusio import read_input, write_output, _get_fillna_dict

from pegasus import tools, misc


import logging
logger = logging.getLogger("pegasus")



def analyze_one_modality(unidata: UnimodalData, output_name: str, is_raw: bool, append_data: UnimodalData, **kwargs) -> None:
    print()
    logger.info(f"Begin to analyze UnimodalData {unidata.get_uid()}.")
    if kwargs["channel_attr"] is not None:
        unidata.obs["Channel"] = unidata.obs[kwargs["channel_attr"]]

    if is_raw:
        # normailize counts and then transform to log space
        tools.log_norm(unidata, kwargs["norm_count"])
        # set group attribute
        if kwargs["batch_correction"] and kwargs["group_attribute"] is not None:
            tools.set_group_attribute(unidata, kwargs["group_attribute"])

    # select highly variable features
    standardize = False # if no select HVF, False
    if kwargs["select_hvf"]:
        if unidata.shape[1] <= kwargs["hvf_ngenes"]:
            logger.warning(f"Number of genes {unidata.shape[1]} is no greater than the target number of highly variable features {kwargs['hvf_ngenes']}. HVF selection is omitted.")
        else:
            standardize = True
            tools.highly_variable_features(
                unidata,
                kwargs["batch_correction"],
                flavor=kwargs["hvf_flavor"],
                n_top=kwargs["hvf_ngenes"],
                n_jobs=kwargs["n_jobs"],
            )
            if kwargs["hvf_flavor"] == "pegasus":
                if kwargs["plot_hvf"] is not None:
                    from pegasus.plotting import hvfplot
                    fig = hvfplot(unidata, show = False)
                    fig.savefig(f"{kwargs['plot_hvf']}.hvf.pdf")

    # batch correction: L/S
    if kwargs["batch_correction"] and kwargs["correction_method"] == "L/S":
        tools.correct_batch(unidata, features="highly_variable_features")

    if kwargs["calc_sigscore"] is not None:
        sig_files = kwargs["calc_sigscore"].split(",")
        for sig_file in sig_files:
            tools.calc_signature_score(unidata, sig_file)

    n_pc = min(kwargs["pca_n"], unidata.shape[0], unidata.shape[1])
    if n_pc < kwargs["pca_n"]:
        logger.warning(f"UnimodalData {unidata.get_uid()} has either dimension ({unidata.shape[0]}, {unidata.shape[1]}) less than the specified number of PCs {kwargs['pca_n']}. Reduce the number of PCs to {n_pc}.")


    if kwargs["batch_correction"] and kwargs["correction_method"] == "scanorama":
        pca_key = tools.run_scanorama(unidata, n_components=n_pc, features="highly_variable_features", standardize=standardize, random_state=kwargs["random_state"])
    else:
        # PCA
        tools.pca(
            unidata,
            n_components=n_pc,
            features="highly_variable_features",
            standardize=standardize,
            robust=kwargs["pca_robust"],
            random_state=kwargs["random_state"],
        )
        pca_key = "pca"

    # batch correction: Harmony
    if kwargs["batch_correction"] and kwargs["correction_method"] == "harmony":
        pca_key = tools.run_harmony(unidata, rep="pca", n_jobs=kwargs["n_jobs"], n_clusters=kwargs["harmony_nclusters"], random_state = kwargs["random_state"])

    # Find K neighbors
    tools.neighbors(
        unidata,
        K=kwargs["K"],
        rep=pca_key,
        n_jobs=kwargs["n_jobs"],
        random_state=kwargs["random_state"],
        full_speed=kwargs["full_speed"],
    )

    # calculate diffmap
    if (
        kwargs["fle"]
        or kwargs["net_fle"]
    ):
        if not kwargs["diffmap"]:
            print("Turn on --diffmap option!")
        kwargs["diffmap"] = True

    if kwargs["diffmap"]:
        tools.diffmap(
            unidata,
            n_components=kwargs["diffmap_ndc"],
            rep=pca_key,
            solver=kwargs["diffmap_solver"],
            random_state=kwargs["random_state"],
            max_t=kwargs["diffmap_maxt"],
        )
        if kwargs["diffmap_to_3d"]:
            tools.reduce_diffmap_to_3d(unidata, random_state=kwargs["random_state"])

    # calculate kBET
    if ("kBET" in kwargs) and kwargs["kBET"]:
        stat_mean, pvalue_mean, accept_rate = tools.calc_kBET(
            unidata,
            kwargs["kBET_batch"],
            rep=pca_key,
            K=kwargs["kBET_K"],
            alpha=kwargs["kBET_alpha"],
            n_jobs=kwargs["n_jobs"],
            random_state=kwargs["random_state"]
        )
        print(
            "kBET stat_mean = {:.2f}, pvalue_mean = {:.4f}, accept_rate = {:.2%}.".format(
                stat_mean, pvalue_mean, accept_rate
            )
        )

    # clustering
    if kwargs["spectral_louvain"]:
        tools.cluster(
            unidata,
            algo="spectral_louvain",
            rep=pca_key,
            resolution=kwargs["spectral_louvain_resolution"],
            rep_kmeans=kwargs["spectral_louvain_basis"],
            n_clusters=kwargs["spectral_louvain_nclusters"],
            n_clusters2=kwargs["spectral_louvain_nclusters2"],
            n_init=kwargs["spectral_louvain_ninit"],
            n_jobs=kwargs["n_jobs"],
            random_state=kwargs["random_state"],
            class_label="spectral_louvain_labels",
        )

    if kwargs["spectral_leiden"]:
        tools.cluster(
            unidata,
            algo="spectral_leiden",
            rep=pca_key,
            resolution=kwargs["spectral_leiden_resolution"],
            rep_kmeans=kwargs["spectral_leiden_basis"],
            n_clusters=kwargs["spectral_leiden_nclusters"],
            n_clusters2=kwargs["spectral_leiden_nclusters2"],
            n_init=kwargs["spectral_leiden_ninit"],
            n_jobs=kwargs["n_jobs"],
            random_state=kwargs["random_state"],
            class_label="spectral_leiden_labels",
        )

    if kwargs["louvain"]:
        tools.cluster(
            unidata,
            algo="louvain",
            rep=pca_key,
            resolution=kwargs["louvain_resolution"],
            random_state=kwargs["random_state"],
            class_label=kwargs["louvain_class_label"],
        )

    if kwargs["leiden"]:
        tools.cluster(
            unidata,
            algo="leiden",
            rep=pca_key,
            resolution=kwargs["leiden_resolution"],
            n_iter=kwargs["leiden_niter"],
            random_state=kwargs["random_state"],
            class_label=kwargs["leiden_class_label"],
        )

    # visualization
    if kwargs["net_tsne"]:
        tools.net_tsne(
            unidata,
            rep=pca_key,
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
            select_frac=kwargs["net_ds_frac"],
            select_K=kwargs["net_ds_K"],
            select_alpha=kwargs["net_ds_alpha"],
            net_alpha=kwargs["net_l2"],
            polish_learning_frac=kwargs["net_tsne_polish_learing_frac"],
            polish_n_iter=kwargs["net_tsne_polish_niter"],
            out_basis=kwargs["net_tsne_basis"],
        )

    if kwargs["net_umap"]:
        tools.net_umap(
            unidata,
            rep=pca_key,
            n_jobs=kwargs["n_jobs"],
            n_neighbors=kwargs["umap_K"],
            min_dist=kwargs["umap_min_dist"],
            spread=kwargs["umap_spread"],
            random_state=kwargs["random_state"],
            select_frac=kwargs["net_ds_frac"],
            select_K=kwargs["net_ds_K"],
            select_alpha=kwargs["net_ds_alpha"],
            full_speed=kwargs["full_speed"],
            net_alpha=kwargs["net_l2"],
            polish_learning_rate=kwargs["net_umap_polish_learing_rate"],
            polish_n_epochs=kwargs["net_umap_polish_nepochs"],
            out_basis=kwargs["net_umap_basis"],
        )

    if kwargs["net_fle"]:
        tools.net_fle(
            unidata,
            output_name,
            n_jobs=kwargs["n_jobs"],
            K=kwargs["fle_K"],
            full_speed=kwargs["full_speed"],
            target_change_per_node=kwargs["fle_target_change_per_node"],
            target_steps=kwargs["fle_target_steps"],
            is3d=False,
            memory=kwargs["fle_memory"],
            random_state=kwargs["random_state"],
            select_frac=kwargs["net_ds_frac"],
            select_K=kwargs["net_ds_K"],
            select_alpha=kwargs["net_ds_alpha"],
            net_alpha=kwargs["net_l2"],
            polish_target_steps=kwargs["net_fle_polish_target_steps"],
            out_basis=kwargs["net_fle_basis"],
        )

    if kwargs["tsne"]:
        tools.tsne(
            unidata,
            rep=pca_key,
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
        )

    if kwargs["fitsne"]:
        tools.fitsne(
            unidata,
            rep=pca_key,
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
        )

    if kwargs["umap"]:
        tools.umap(
            unidata,
            rep=pca_key,
            n_neighbors=kwargs["umap_K"],
            min_dist=kwargs["umap_min_dist"],
            spread=kwargs["umap_spread"],
            random_state=kwargs["random_state"],
        )

    if kwargs["fle"]:
        tools.fle(
            unidata,
            output_name,
            n_jobs=kwargs["n_jobs"],
            K=kwargs["fle_K"],
            full_speed=kwargs["full_speed"],
            target_change_per_node=kwargs["fle_target_change_per_node"],
            target_steps=kwargs["fle_target_steps"],
            is3d=False,
            memory=kwargs["fle_memory"],
            random_state=kwargs["random_state"],
        )

    # calculate diffusion-based pseudotime from roots
    if len(kwargs["pseudotime"]) > 0:
        tools.calc_pseudotime(unidata, kwargs["pseudotime"])


    genome = unidata.uns["genome"]

    if append_data is not None:
        locs = unidata.obs_names.get_indexer(append_data.obs_names)
        idx = locs >= 0
        locs = locs[idx]
        Y = append_data.X[idx, :].tocoo(copy = False)
        Z = coo_matrix((Y.data, (locs[Y.row], Y.col)), shape = (unidata.shape[0], append_data.shape[1])).tocsr()

        idy = Z.getnnz(axis = 0) > 0
        n_nonzero = idy.sum()
        if n_nonzero > 0:
            if n_nonzero < append_data.shape[1]:
                Z = Z[:, idy]
                append_df = append_data.feature_metadata.loc[idy, :]
            else:
                append_df = append_data.feature_metadata

            rawX = hstack([unidata.get_matrix("raw.X"), Z], format = "csr")

            Zt = Z.astype(np.float32)
            Zt.data *= np.repeat(unidata.obs["scale"].values, np.diff(Zt.indptr))
            Zt.data = np.log1p(Zt.data)

            X = hstack([unidata.get_matrix("X"), Zt], format = "csr")

            new_genome = unidata.get_genome() + "_and_" + append_data.get_genome()

            feature_metadata = pd.concat([unidata.feature_metadata, append_df], axis = 0)
            feature_metadata.reset_index(inplace = True)
            feature_metadata.fillna(value = _get_fillna_dict(unidata.feature_metadata), inplace = True)

            unidata = UnimodalData(unidata.barcode_metadata, feature_metadata, {"X": X, "raw.X": rawX}, unidata.uns.mapping, unidata.obsm.mapping, unidata.varm.mapping) # uns.mapping, obsm.mapping and varm.mapping are passed by reference
            unidata.uns["genome"] = new_genome


    if kwargs["output_h5ad"]:
        adata = unidata.to_anndata()
        adata.uns["scale.data"] = adata.uns.pop("fmat_highly_variable_features")  # assign by reference
        adata.uns["scale.data.rownames"] = unidata.var_names[unidata.var["highly_variable_features"]].values
        adata.write(f"{output_name}.h5ad", compression="gzip")
        del adata

    # write out results
    if kwargs["output_loom"]:
        write_output(unidata, f"{output_name}.loom")


    # Change genome name back if append_data is True
    if unidata.uns["genome"] != genome:
        unidata.uns["genome"] = genome
    # Eliminate objects starting with fmat_ from uns
    for key in list(unidata.uns):
        if key.startswith("fmat_"):
            unidata.uns.pop(key)



def run_pipeline(input_file: str, output_name: str, **kwargs):
    is_raw = not kwargs["processed"]

    black_list = set()
    if kwargs["black_list"] is not None:
        black_list = set(kwargs["black_list"].split(","))

    # load input data
    data = read_input(input_file, black_list = black_list)

    # process focus_list
    focus_list = kwargs["focus"]
    if len(focus_list) == 0:
        focus_list = [data.current_key()]


    append_data = None
    if kwargs["append"] is not None:
        append_data = data.get_data(kwargs["append"])

    logger.info("Inputs are loaded.")

    if is_raw and not kwargs["subcluster"]:
        # filter out low quality cells/genes
        tools._run_filter_data(
            data,
            focus_list=focus_list,
            output_filt=kwargs["output_filt"],
            plot_filt=kwargs["plot_filt"],
            plot_filt_figsize=kwargs["plot_filt_figsize"],
            min_genes_before_filt=kwargs["min_genes_before_filt"],
            select_singlets=kwargs["select_singlets"],
            remap_string=kwargs["remap_singlets"],
            subset_string=kwargs["subset_singlets"],
            min_genes=kwargs["min_genes"],
            max_genes=kwargs["max_genes"],
            min_umis=kwargs["min_umis"],
            max_umis=kwargs["max_umis"],
            mito_prefix=kwargs["mito_prefix"],
            percent_mito=kwargs["percent_mito"],
            percent_cells=kwargs["percent_cells"],
        )


    for key in focus_list:
        unidata = data.get_data(key)
        analyze_one_modality(unidata, f"{output_name}.{unidata.get_uid()}", is_raw, append_data, **kwargs)

    print()

    # if kwargs["subcluster"]:
    #     unidata = tools.get_anndata_for_subclustering(adata, kwargs["subset_selections"])
    #     is_raw = True  # get submat and then set is_raw to True

    # write out results

    write_output(data, f"{output_name}.zarr.zip")

    print("Results are written.")
