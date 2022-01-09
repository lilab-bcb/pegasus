import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, coo_matrix, hstack

from pegasusio import UnimodalData, MultimodalData
from pegasusio import read_input, write_output, _fillna

from pegasus import tools, misc


import logging
logger = logging.getLogger("pegasus")



def analyze_one_modality(unidata: UnimodalData, output_name: str, is_raw: bool, append_data: UnimodalData, **kwargs) -> None:
    print()
    logger.info(f"Begin to analyze UnimodalData {unidata.get_uid()}.")

    if is_raw:
        # normailize counts and then transform to log space
        tools.log_norm(unidata, kwargs["norm_count"])

        # select highly variable features
        standardize = False # if no select HVF, False
        if kwargs["select_hvf"]:
            if unidata.shape[1] <= kwargs["hvf_ngenes"]:
                logger.warning(f"Number of genes {unidata.shape[1]} is no greater than the target number of highly variable features {kwargs['hvf_ngenes']}. HVF selection is omitted.")
            else:
                standardize = True
                tools.highly_variable_features(
                    unidata,
                    kwargs["batch_attr"] if kwargs["batch_correction"] else None,
                    flavor=kwargs["hvf_flavor"],
                    n_top=kwargs["hvf_ngenes"],
                    n_jobs=kwargs["n_jobs"],
                )
                if kwargs["hvf_flavor"] == "pegasus":
                    if kwargs["plot_hvf"] is not None:
                        from pegasus.plotting import hvfplot
                        fig = hvfplot(unidata, return_fig = True)
                        fig.savefig(f"{kwargs['plot_hvf']}.hvf.pdf")


        n_pc = min(kwargs["pca_n"], unidata.shape[0], unidata.shape[1])
        if n_pc < kwargs["pca_n"]:
            logger.warning(f"UnimodalData {unidata.get_uid()} has either dimension ({unidata.shape[0]}, {unidata.shape[1]}) less than the specified number of PCs {kwargs['pca_n']}. Reduce the number of PCs to {n_pc}.")

        # Run PCA irrespect of which batch correction method would apply
        tools.pca(
            unidata,
            n_components=n_pc,
            features="highly_variable_features",
            standardize=standardize,
            n_jobs=kwargs["n_jobs"],
            random_state=kwargs["random_state"],
        )
        dim_key = "pca"



        if kwargs["nmf"] or (kwargs["batch_correction"] and kwargs["correction_method"] == "inmf"):
            n_nmf = min(kwargs["nmf_n"], unidata.shape[0], unidata.shape[1])
            if n_nmf < kwargs["nmf_n"]:
                logger.warning(f"UnimodalData {unidata.get_uid()} has either dimension ({unidata.shape[0]}, {unidata.shape[1]}) less than the specified number of NMF components {kwargs['nmf_n']}. Reduce the number of NMF components to {n_nmf}.")

        if kwargs["nmf"]:
            if kwargs["batch_correction"] and kwargs["correction_method"] == "inmf":
                logger.warning("NMF is skipped because integrative NMF is run instead.")
            else:
                tools.nmf(
                    unidata,
                    n_components=n_nmf,
                    features="highly_variable_features",
                    n_jobs=kwargs["n_jobs"],
                    random_state=kwargs["random_state"],
                )


        if kwargs["batch_correction"]:
            if kwargs["correction_method"] == "harmony":
                dim_key = tools.run_harmony(unidata, batch=kwargs["batch_attr"], rep="pca", n_jobs=kwargs["n_jobs"], n_clusters=kwargs["harmony_nclusters"], random_state = kwargs["random_state"])
            elif kwargs["correction_method"] == "inmf":
                dim_key = tools.integrative_nmf(unidata, batch=kwargs["batch_attr"], n_components=n_nmf, features="highly_variable_features", lam=kwargs["inmf_lambda"], n_jobs=kwargs["n_jobs"], random_state = kwargs["random_state"])
            elif kwargs["correction_method"] == "scanorama":
                dim_key = tools.run_scanorama(unidata, batch=kwargs["batch_attr"], n_components=n_pc, features="highly_variable_features", standardize=standardize, random_state=kwargs["random_state"])
            else:
                raise ValueError(f"Unknown batch correction method {kwargs['correction_method']}!")


        # Find K neighbors
        tools.neighbors(
            unidata,
            K=kwargs["K"],
            rep=dim_key,
            n_jobs=kwargs["n_jobs"],
            random_state=kwargs["random_state"],
            full_speed=kwargs["full_speed"],
        )


    if kwargs["calc_sigscore"] is not None:
        sig_files = kwargs["calc_sigscore"].split(",")
        for sig_file in sig_files:
            tools.calc_signature_score(unidata, sig_file)

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
            rep=dim_key,
            solver=kwargs["diffmap_solver"],
            max_t=kwargs["diffmap_maxt"],
            n_jobs=kwargs["n_jobs"],
            random_state=kwargs["random_state"],
        )

    # calculate kBET
    if ("kBET" in kwargs) and kwargs["kBET"]:
        stat_mean, pvalue_mean, accept_rate = tools.calc_kBET(
            unidata,
            kwargs["kBET_batch"],
            rep=dim_key,
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
            rep=dim_key,
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
            rep=dim_key,
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
            rep=dim_key,
            resolution=kwargs["louvain_resolution"],
            random_state=kwargs["random_state"],
            class_label=kwargs["louvain_class_label"],
        )

    if kwargs["leiden"]:
        tools.cluster(
            unidata,
            algo="leiden",
            rep=dim_key,
            resolution=kwargs["leiden_resolution"],
            n_iter=kwargs["leiden_niter"],
            random_state=kwargs["random_state"],
            class_label=kwargs["leiden_class_label"],
        )

    # visualization
    if kwargs["net_umap"]:
        tools.net_umap(
            unidata,
            rep=dim_key,
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
            rep=dim_key,
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
            initialization=kwargs["tsne_init"],
        )

    if kwargs["umap"]:
        tools.umap(
            unidata,
            rep=dim_key,
            n_neighbors=kwargs["umap_K"],
            min_dist=kwargs["umap_min_dist"],
            spread=kwargs["umap_spread"],
            n_jobs=kwargs["n_jobs"],
            full_speed=kwargs["full_speed"],
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

    if kwargs["infer_doublets"]:
        channel_attr = "Channel"
        if (channel_attr not in unidata.obs) or (unidata.obs["Channel"].cat.categories.size == 1):
            channel_attr = None
        clust_attr = kwargs["dbl_cluster_attr"]
        if (clust_attr is None) or (clust_attr not in unidata.obs):
            clust_attr = None
            for value in ["leiden_labels", "louvain_labels", "spectral_leiden_labels", "spectral_louvain_labels"]:
                if value in unidata.obs:
                    clust_attr = value
                    break

        if channel_attr is not None:
            logger.info(f"For doublet inference, channel_attr={channel_attr}.")
        if clust_attr is not None:
            logger.info(f"For doublet inference, clust_attr={clust_attr}.")

        tools.infer_doublets(unidata, channel_attr = channel_attr, clust_attr = clust_attr, expected_doublet_rate = kwargs["expected_doublet_rate"], n_jobs = kwargs["n_jobs"], random_state = kwargs["random_state"], plot_hist = output_name)

        dbl_clusts = None
        if clust_attr is not None:
            clusts = []
            for idx, row in unidata.uns["pred_dbl_cluster"].iterrows():
                if row["percentage"] >= 50.0:
                    logger.info(f"Cluster {row['cluster']} (percentage={row['percentage']:.2f}%, q-value={row['qval']:.6g}) is identified as a doublet cluster.")
                    clusts.append(row["cluster"])
            if len(clusts) > 0:
                dbl_clusts = f"{clust_attr}:{','.join(clusts)}"

        tools.mark_doublets(unidata, dbl_clusts = dbl_clusts)


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

            if kwargs["citeseq"]:
                append_df = append_df.copy()
                append_df.index = append_df.index.map(lambda x: f"Ab-{x}")

            rawX = hstack([unidata.get_matrix("raw.X"), Z], format = "csr")

            Zt = Z.astype(np.float32)
            if not kwargs["citeseq"]:
                Zt.data *= np.repeat(unidata.obs["scale"].values, np.diff(Zt.indptr))
                Zt.data = np.log1p(Zt.data)
            else:
                Zt.data = np.arcsinh(Zt.data / 5.0, dtype = np.float32)

            X = hstack([unidata.get_matrix("X"), Zt], format = "csr")

            new_genome = unidata.get_genome()
            if new_genome != append_data.get_genome():
                new_genome = f"{new_genome}_and_{append_data.get_genome()}"

            feature_metadata = pd.concat([unidata.feature_metadata, append_df], axis = 0)
            feature_metadata.reset_index(inplace = True)
            _fillna(feature_metadata)
            unidata = UnimodalData(unidata.barcode_metadata, feature_metadata, {"X": X, "raw.X": rawX}, unidata.uns.mapping, unidata.obsm.mapping, unidata.varm.mapping) # uns.mapping, obsm.mapping and varm.mapping are passed by reference
            unidata.uns["genome"] = new_genome

            if kwargs["citeseq"] and kwargs["citeseq_umap"]:
                umap_index = append_df.index.difference([f"Ab-{x}" for x in kwargs["citeseq_umap_exclude"]])
                unidata.obsm["X_citeseq"] = unidata.X[:, unidata.var_names.isin(umap_index)].toarray()
                tools.umap(
                    unidata,
                    rep="citeseq",
                    n_neighbors=kwargs["umap_K"],
                    min_dist=kwargs["umap_min_dist"],
                    spread=kwargs["umap_spread"],
                    n_jobs=kwargs["n_jobs"],
                    full_speed=kwargs["full_speed"],
                    random_state=kwargs["random_state"],
                    out_basis="citeseq_umap",
                )

    if kwargs["output_h5ad"]:
        import time
        start_time = time.perf_counter()
        adata = unidata.to_anndata()
        if "_tmp_fmat_highly_variable_features" in adata.uns:
            adata.uns["scale.data"] = adata.uns.pop("_tmp_fmat_highly_variable_features")  # assign by reference
            adata.uns["scale.data.rownames"] = unidata.var_names[unidata.var["highly_variable_features"]==True].values
        adata.write(f"{output_name}.h5ad", compression="gzip")
        del adata
        end_time = time.perf_counter()
        logger.info(f"H5AD file {output_name}.h5ad is written. Time spent = {end_time - start_time:.2f}s.")

    # write out results
    if kwargs["output_loom"]:
        write_output(unidata, f"{output_name}.loom")


    # Change genome name back if append_data is True
    if unidata.uns["genome"] != genome:
        unidata.uns["genome"] = genome
    # Eliminate objects starting with _tmp from uns
    unidata.uns.pop("_tmp_fmat_highly_variable_features", None)



def run_pipeline(input_file: str, output_name: str, **kwargs):
    is_raw = not kwargs["processed"]

    black_list = set()
    if kwargs["black_list"] is not None:
        black_list = set(kwargs["black_list"].split(","))

    genome = None
    if kwargs["genome"] is not None:
        genome = kwargs["genome"]

    # load input data
    data = read_input(input_file, black_list=black_list, genome=genome)

    if kwargs["citeseq"]:
        # temporary solution for CITE-Seq data
        keys = list(data.data)
        assert len(keys) == 2
        if keys[0].endswith("-citeseq"):
            keys = [keys[1], keys[0]]
        assert keys[0].endswith("-rna") and keys[1].endswith("-citeseq")
        focus_list = [keys[0]]
        append_data = data.get_data(keys[1])

        kwargs["output_h5ad"] = True
    else:
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
