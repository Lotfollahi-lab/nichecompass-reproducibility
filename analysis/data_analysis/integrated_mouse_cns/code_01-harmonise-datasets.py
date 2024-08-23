import anndata as ad
import os
import pandas as pd
from typing import Literal

# identify shared genes

merfish_anndata = ad.read_h5ad("data/merfish_mouse_brain/merfish_mouse_brain_batch1.h5ad")
merfish_genes = [x.capitalize() for x in merfish_anndata.var_names.tolist()]

starmap_anndata = ad.read_h5ad("data/starmap_plus_mouse_cns/starmap_plus_mouse_cns_batch1.h5ad")
starmap_genes = [x.capitalize() for x in starmap_anndata.var_names.tolist()]

shared_genes = list(set(merfish_genes).intersection(set(starmap_genes)))


# load and standardize the anndata objects

def load():

    def _standardize_anndata(anndata, genes, dataset: Literal["merfish", "starmap"]):

        dataset_section_column = {
            "merfish": "brain_section_label",
            "starmap": "batch"}

        raw_counts = anndata.layers["counts"]
        var_names = [x.capitalize() for x in anndata.var_names]
        obs_names = anndata.obs_names
        obs = pd.DataFrame.from_dict({
            "section": anndata.obs[dataset_section_column[dataset]],
        })
        obs["dataset"] = dataset
        spatial = anndata.obsm["spatial"]

        standardized_anndata = ad.AnnData(
            raw_counts,
            obs=obs,
            obsm={"spatial": spatial},
        )
        standardized_anndata.var_names = var_names
        standardized_anndata.obs_names = obs_names

        standardized_anndata = standardized_anndata[:, genes]

        section = anndata.obs[dataset_section_column[dataset]].tolist()[0]
        return dataset, section, standardized_anndata

    datasets = []

    slide_files = os.listdir("data/merfish_mouse_brain/")
    for i, slide in enumerate(os.listdir("data/merfish_mouse_brain/")):
        print(f"Processing merfish slide {i + 1} of {len(slide_files)}: {slide}")
        anndata = ad.read_h5ad(f"data/merfish_mouse_brain/{slide}")
        dataset, section, standardized_anndata = _standardize_anndata(anndata, shared_genes, "merfish")
        datasets.append((dataset, section, standardized_anndata))

    slide_files = os.listdir("data/starmap_plus_mouse_cns/")
    for i, slide in enumerate(slide_files):
        print(f"Processing starmap slide {i + 1} of {len(slide_files)}: {slide}")
        anndata = ad.read_h5ad(f"data/starmap_plus_mouse_cns/{slide}")
        dataset, section, standardized_anndata = _standardize_anndata(anndata, shared_genes, "starmap")
        datasets.append((dataset, section, standardized_anndata))

    return datasets


datasets = load()
os.makedirs("nichecompass-mouse-cns-integration", exist_ok=True)
for dataset in datasets:
    dataset, section, standardized_anndata = dataset
    file_name = f"{dataset}__{section}.h5ad"
    local_file = os.path.join("nichecompass-mouse-cns-integration", file_name)
    standardized_anndata.write(local_file)