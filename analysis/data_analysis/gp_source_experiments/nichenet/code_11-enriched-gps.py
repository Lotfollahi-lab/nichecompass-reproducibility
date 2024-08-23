import os
import anndata as ad
import pickle

from nichecompass.models import NicheCompass


def enrich_gps(model_directory, anndata):

    model = NicheCompass.load(dir_path=model_directory, adata=anndata)

    anndata_with_enriched_gps = anndata.copy()

    enriched_gps = model.run_differential_gp_tests(
        cat_key="nichecompass_latent_cluster",
        selected_cats=None,
        comparison_cats="rest",
        log_bayes_factor_thresh=2.3,
        adata=anndata_with_enriched_gps)

    return enriched_gps, anndata_with_enriched_gps


trained_model_artifact_directory = "trained-model:nichenet"
trained_model_umaps_directory = "trained-model-umaps:nichenet"

anndata = ad.read_h5ad(os.path.join(trained_model_umaps_directory, "anndata_umap_with_clusters.h5ad"))

enriched_gps, anndata_with_enriched_gps = enrich_gps(trained_model_artifact_directory, anndata)

os.makedirs("enriched-gps:nichenet", exist_ok=True)

file_name = f"enriched_gps.pkl"
local_file = os.path.join("enriched-gps:nichenet", file_name)
with open(local_file, "wb") as file:
    pickle.dump(enriched_gps, file)

file_name = "anndata_with_enriched_gps.h5ad"
local_file = os.path.join("enriched-gps:nichenet", file_name)
anndata_with_enriched_gps.write_h5ad(local_file)

