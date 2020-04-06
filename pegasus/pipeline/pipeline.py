import anndata
import numpy as np
from scipy.sparse import csr_matrix, hstack

from pegasus import io, tools, cite_seq, misc


def run_pipeline(input_file, output_name, **kwargs):
    is_raw = not kwargs["processed"]

    if "seurat_compatible" not in kwargs:
        kwargs["seurat_compatible"] = False

    # load input data
    adata = io.read_input(
        input_file,
        genome=kwargs["genome"],
        concat_matrices=False if kwargs["cite_seq"] else True,
        h5ad_mode=("a" if (is_raw or kwargs["subcluster"]) else "r+"),
        select_singlets=kwargs["select_singlets"],
        channel_attr=kwargs["channel_attr"],
        black_list=(
            kwargs["black_list"].split(",") if kwargs["black_list"] is not None else []
        ),
    )

    if kwargs["cite_seq"]:
        data_list = adata
        assert len(data_list) == 2
        adata = cdata = None
        for i in range(len(data_list)):
            if data_list[i].uns["genome"].startswith("CITE_Seq"):
                cdata = data_list[i]
            else:
                adata = data_list[i]
        assert adata is not None and cdata is not None

    if is_raw:
        values = adata.X.getnnz(axis=1)
        if values.min() == 0:  # 10x raw data
            adata._inplace_subset_obs(values >= kwargs["min_genes_on_raw"])
        if kwargs['remap_singlets'] is not None:
            misc.remap_singlets(adata, kwargs['remap_singlets'])
        if kwargs['subset_singlets'] is not None:
            misc.subset_singlets(adata, kwargs['subset_singlets'])

    print("Inputs are loaded.")

    if kwargs["seurat_compatible"]:
        assert is_raw and kwargs["select_hvf"]

    if kwargs["subcluster"]:
        adata = tools.get_anndata_for_subclustering(adata, kwargs["subset_selections"])
        is_raw = True  # get submat and then set is_raw to True

    if is_raw:
        if not kwargs["subcluster"]:
            # filter out low quality cells/genes
            tools.run_filter_data(
                adata,
                output_filt=kwargs["output_filt"],
                plot_filt=kwargs["plot_filt"],
                plot_filt_figsize=kwargs["plot_filt_figsize"],
                mito_prefix=kwargs["mito_prefix"],
                min_genes=kwargs["min_genes"],
                max_genes=kwargs["max_genes"],
                min_umis=kwargs["min_umis"],
                max_umis=kwargs["max_umis"],
                percent_mito=kwargs["percent_mito"],
                percent_cells=kwargs["percent_cells"],
            )

            if kwargs["seurat_compatible"]:
                raw_data = adata.copy()  # raw as count

            # normailize counts and then transform to log space
            tools.log_norm(adata, kwargs["norm_count"])

            # set group attribute
            if kwargs["batch_correction"] and kwargs["group_attribute"] is not None:
                tools.set_group_attribute(adata, kwargs["group_attribute"])

        # select highly variable features
        if kwargs["select_hvf"]:
            tools.highly_variable_features(
                adata,
                kwargs["batch_correction"],
                flavor=kwargs["hvf_flavor"],
                n_top=kwargs["hvf_ngenes"],
                n_jobs=kwargs["n_jobs"],
            )
            if kwargs["hvf_flavor"] == "pegasus":
                if kwargs["plot_hvf"] is not None:
                    from pegasus.plotting import plot_hvf

                    robust_idx = adata.var["robust"].values
                    plot_hvf(
                        adata.var.loc[robust_idx, "mean"],
                        adata.var.loc[robust_idx, "var"],
                        adata.var.loc[robust_idx, "hvf_loess"],
                        adata.var.loc[robust_idx, "highly_variable_features"],
                        kwargs["plot_hvf"] + ".hvf.pdf",
                    )

        # batch correction: L/S
        if kwargs["batch_correction"] and kwargs["correction_method"] == "L/S":
            tools.correct_batch(adata, features="highly_variable_features")

        # PCA
        tools.pca(
            adata,
            n_components=kwargs["pca_n"],
            features="highly_variable_features",
            robust=kwargs["pca_robust"],
            random_state=kwargs["random_state"],
        )

        pca_key = "pca"
        # batch correction: Harmony
        if kwargs["batch_correction"] and kwargs["correction_method"] == "harmony":
            pca_key = tools.run_harmony(adata, rep="pca", n_jobs=kwargs["n_jobs"], n_clusters=kwargs["harmony_nclusters"], random_state = kwargs["random_state"])

        # Find K neighbors
        tools.neighbors(
            adata,
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
                adata,
                n_components=kwargs["diffmap_ndc"],
                rep=pca_key,
                solver=kwargs["diffmap_solver"],
                random_state=kwargs["random_state"],
                max_t=kwargs["diffmap_maxt"],
            )
            if kwargs["diffmap_to_3d"]:
                tools.reduce_diffmap_to_3d(adata, random_state=kwargs["random_state"])

    # calculate kBET
    if ("kBET" in kwargs) and kwargs["kBET"]:
        stat_mean, pvalue_mean, accept_rate = tools.calc_kBET(
            adata,
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
            adata,
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
            adata,
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
            adata,
            algo="louvain",
            rep=pca_key,
            resolution=kwargs["louvain_resolution"],
            random_state=kwargs["random_state"],
            class_label=kwargs["louvain_class_label"],
        )

    if kwargs["leiden"]:
        tools.cluster(
            adata,
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
            adata,
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
            adata,
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
            adata,
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
            adata,
            rep=pca_key,
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
        )

    if kwargs["fitsne"]:
        tools.fitsne(
            adata,
            rep=pca_key,
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
        )

    if kwargs["umap"]:
        tools.umap(
            adata,
            rep=pca_key,
            n_neighbors=kwargs["umap_K"],
            min_dist=kwargs["umap_min_dist"],
            spread=kwargs["umap_spread"],
            random_state=kwargs["random_state"],
        )

    if kwargs["fle"]:
        tools.fle(
            adata,
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
        tools.calc_pseudotime(adata, kwargs["pseudotime"])

    # merge cite-seq data and run t-SNE
    if kwargs["cite_seq"]:
        adt_matrix = np.zeros((adata.shape[0], cdata.shape[1]), dtype="float32")
        idx = adata.obs_names.isin(cdata.obs_names)
        adt_matrix[idx, :] = cdata[adata.obs_names[idx],].X.toarray()
        if abs(100.0 - kwargs["cite_seq_capping"]) > 1e-4:
            cite_seq.capping(adt_matrix, kwargs["cite_seq_capping"])

        var_names = np.concatenate(
            [adata.var_names, ["AD-" + x for x in cdata.var_names]]
        )

        new_data = anndata.AnnData(
            X=hstack([adata.X, csr_matrix(adt_matrix)], format="csr"),
            obs=adata.obs,
            obsm=adata.obsm,
            uns=adata.uns,
            var={
                "var_names": var_names,
                "gene_ids": var_names,
                "n_cells": np.concatenate(
                    [adata.var["n_cells"].values, [0] * cdata.shape[1]]
                ),
                "percent_cells": np.concatenate(
                    [adata.var["percent_cells"].values, [0.0] * cdata.shape[1]]
                ),
                "robust": np.concatenate(
                    [adata.var["robust"].values, [False] * cdata.shape[1]]
                ),
                "highly_variable_features": np.concatenate(
                    [
                        adata.var["highly_variable_features"].values,
                        [False] * cdata.shape[1],
                    ]
                ),
            },
        )
        new_data.obsm["X_CITE-Seq"] = adt_matrix
        adata = new_data
        print("ADT count matrix is attached.")

        tools.fitsne(
            adata,
            rep="CITE-Seq",
            n_jobs=kwargs["n_jobs"],
            perplexity=kwargs["tsne_perplexity"],
            random_state=kwargs["random_state"],
            out_basis="citeseq_fitsne",
        )
        print("Antibody embedding is done.")

    if kwargs["seurat_compatible"]:
        seurat_data = adata.copy()
        seurat_data.raw = raw_data
        seurat_data.uns["scale.data"] = adata.uns["fmat_highly_variable_features"]  # assign by reference
        seurat_data.uns["scale.data.rownames"] = adata.var_names[
            adata.var["highly_variable_features"]
        ].values
        io.write_output(seurat_data, output_name + ".seurat.h5ad")

    # write out results
    io.write_output(adata, output_name + ".h5ad")

    if kwargs["output_loom"]:
        io.write_output(adata, output_name + ".loom")

    print("Results are written.")
