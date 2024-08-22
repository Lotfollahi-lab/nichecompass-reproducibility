import scanpy as sc
import os
import anndata as ad


config = {
   "resolution": 0.4
}


def cluster(anndata, leiden_resolution):

    sc.settings.n_jobs = -1

    anndata_with_clusters = sc.tl.leiden(
        adata=anndata,
        resolution=0.4,
        neighbors_key="nichecompass_latent",
        key_added="nichecompass_latent_cluster",
        copy=True)

    return anndata_with_clusters


def generate_umaps(anndata_with_clusters):

    umap_by_dataset_fig =  sc.pl.umap(anndata_with_clusters, color="dataset", size=2, return_fig=True)
    umap_by_section_fig = sc.pl.umap(anndata_with_clusters, color="section", size=2, return_fig=True)
    umap_by_niche_fig = sc.pl.umap(anndata_with_clusters, color="nichecompass_latent_cluster", size=2, return_fig=True)

    return {"by_dataset": umap_by_dataset_fig, "by_section": umap_by_section_fig, "by_niche": umap_by_niche_fig}


trained_model_artifact_directory = "trained-model"
trained_model_anndata = ad.read(os.path.join(trained_model_artifact_directory, "adata.h5ad"))
anndata_with_clusters = cluster(trained_model_anndata, config["resolution"])
figures = generate_umaps(anndata_with_clusters)

os.makedirs("trained-model-umaps", exist_ok=True)

for figure_name, figure in figures.items():
    file_name = f"umap_{figure_name}.png"
    local_file = os.path.join("trained-model-umaps", file_name)
    figure.savefig(local_file, bbox_inches='tight', dpi=200)

file_name = f"anndata_umap_with_clusters.h5ad"
local_file = os.path.join("trained-model-umaps", file_name)
anndata_with_clusters.write_h5ad(local_file)
