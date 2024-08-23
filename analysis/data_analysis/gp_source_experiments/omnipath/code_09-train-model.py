import scanpy as sc
import os
import anndata as ad

from nichecompass.models import NicheCompass


config = {
    "cat_covariates_keys": ["sample"],
    "cat_covariates_embeds_injection": ["gene_expr_decoder"],
    "cat_covariates_embeds_nums": [3],
    "cat_covariates_no_edges": [True],
    "conv_layer_encoder": "gatv2conv",
    "active_gp_thresh_ratio": 0.01,
    "n_epochs": 100,
    "n_epochs_all_gps": 25,
    "lr": 0.001,
    "lambda_edge_recon": 5000000.,
    "lambda_gene_expr_recon": 3000.,
    "lambda_l1_masked": 0.,
    "lambda_l1_addon": 0.,
    "edge_batch_size": 256,
    "n_sampled_neighbors": 4,
    "n_addon_gp": 100,
}


def train(
        adata,
        output_directory,
        cat_covariates_keys,
        cat_covariates_embeds_injection,
        cat_covariates_embeds_nums,
        cat_covariates_no_edges,
        conv_layer_encoder,
        active_gp_thresh_ratio,
        n_epochs,
        n_epochs_all_gps,
        lr,
        lambda_edge_recon,
        lambda_gene_expr_recon,
        lambda_l1_masked,
        lambda_l1_addon,
        edge_batch_size,
        n_sampled_neighbors,
        n_addon_gp
):

    adj_key = "spatial_connectivities"
    gp_names_key = "nichecompass_gp_names"
    active_gp_names_key = "nichecompass_active_gp_names"
    gp_targets_mask_key = "nichecompass_gp_targets"
    gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
    gp_sources_mask_key = "nichecompass_gp_sources"
    gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
    latent_key = "nichecompass_latent"

    use_cuda_if_available = True

    # Initialize model
    model = NicheCompass(adata,
                         counts_key=None,
                         adj_key=adj_key,
                         cat_covariates_embeds_injection=cat_covariates_embeds_injection,
                         cat_covariates_keys=cat_covariates_keys,
                         cat_covariates_no_edges=cat_covariates_no_edges,
                         cat_covariates_embeds_nums=cat_covariates_embeds_nums,
                         gp_names_key=gp_names_key,
                         active_gp_names_key=active_gp_names_key,
                         gp_targets_mask_key=gp_targets_mask_key,
                         gp_targets_categories_mask_key=gp_targets_categories_mask_key,
                         gp_sources_mask_key=gp_sources_mask_key,
                         gp_sources_categories_mask_key=gp_sources_categories_mask_key,
                         latent_key=latent_key,
                         conv_layer_encoder=conv_layer_encoder,
                         active_gp_thresh_ratio=active_gp_thresh_ratio,
                         n_addon_gp=n_addon_gp)

    # Train model
    model.train(n_epochs=n_epochs,
                n_epochs_all_gps=n_epochs_all_gps,
                lr=lr,
                lambda_edge_recon=lambda_edge_recon,
                lambda_gene_expr_recon=lambda_gene_expr_recon,
                lambda_l1_masked=lambda_l1_masked,
                lambda_l1_addon=lambda_l1_addon,
                edge_batch_size=edge_batch_size,
                n_sampled_neighbors=n_sampled_neighbors,
                use_cuda_if_available=use_cuda_if_available,
                verbose=False)

    # Compute latent neighbor graph
    sc.pp.neighbors(model.adata,
                    use_rep=latent_key,
                    key_added=latent_key)

    # Compute UMAP embedding
    sc.tl.umap(model.adata,
               neighbors_key=latent_key)

    # Save trained model
    model.save(dir_path=output_directory,
               overwrite=True,
               save_adata=True,
               adata_file_name="adata.h5ad")


dataset_directory = "dataset-with-gene-programs:omnipath"

dataset_anndata = ad.read(os.path.join(dataset_directory, "dataset_with_gene_programs.h5ad"))

os.makedirs("trained-model:omnipath", exist_ok=True)

train(
    adata=dataset_anndata,
    output_directory="trained-model:omnipath",
    cat_covariates_keys=config["cat_covariates_keys"],
    cat_covariates_embeds_injection=config["cat_covariates_embeds_injection"],
    cat_covariates_embeds_nums=config["cat_covariates_embeds_nums"],
    cat_covariates_no_edges=config["cat_covariates_no_edges"],
    conv_layer_encoder=config["conv_layer_encoder"],
    active_gp_thresh_ratio=config["active_gp_thresh_ratio"],
    n_epochs=config["n_epochs"],
    n_epochs_all_gps=config["n_epochs_all_gps"],
    lr=config["lr"],
    lambda_edge_recon=config["lambda_edge_recon"],
    lambda_gene_expr_recon=config["lambda_gene_expr_recon"],
    lambda_l1_masked=config["lambda_l1_masked"],
    lambda_l1_addon=config["lambda_l1_addon"],
    edge_batch_size=config["edge_batch_size"],
    n_sampled_neighbors=config["n_sampled_neighbors"],
    n_addon_gp=config["n_addon_gp"]
)
