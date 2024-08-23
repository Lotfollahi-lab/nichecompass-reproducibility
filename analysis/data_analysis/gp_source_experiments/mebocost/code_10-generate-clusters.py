import os
import anndata as ad
import rapids_singlecell as rsc


config = {
   "resolution": 0.36
}


def cluster(anndata, leiden_resolution):

    print("moving anndata to GPU")
    rsc.get.anndata_to_GPU(anndata)

    print("performing clustering")
    anndata_with_clusters = rsc.tl.leiden(
        adata=anndata,
        resolution=leiden_resolution,
        neighbors_key="nichecompass_latent",
        key_added="nichecompass_latent_cluster",
        copy=True)

    print("moving result anndata to CPU")
    rsc.get.anndata_to_CPU(anndata_with_clusters)

    print("result moved to CPU")
    print(anndata_with_clusters)
    return anndata_with_clusters


trained_model_artifact_directory = "trained-model:mebocost"

trained_model_anndata = ad.read(os.path.join(trained_model_artifact_directory, "adata.h5ad"))
anndata_with_clusters = cluster(trained_model_anndata, config["resolution"])

os.makedirs("trained-model-umaps:mebocost", exist_ok=True)

file_name = f"anndata_umap_with_clusters.h5ad"
local_file = os.path.join("trained-model-umaps:mebocost", file_name)
anndata_with_clusters.write_h5ad(local_file)
