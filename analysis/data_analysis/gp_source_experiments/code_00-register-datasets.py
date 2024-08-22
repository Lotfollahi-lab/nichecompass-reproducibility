import anndata as ad
import os
import pandas as pd

# load and standardize the anndata objects

def load():

    def standardize_anndata(dataset_name, excluded_genes, anndata):

        if excluded_genes is None:
            excluded_genes = []

        gene_names = anndata.var_names.tolist()
        gene_names_filtered = list(set(gene_names).difference(excluded_genes))

        anndata = anndata[:, gene_names_filtered]

        dataset_section_column = "batch"

        raw_counts = anndata.layers["counts"]
        var_names = [x.capitalize() for x in anndata.var_names]
        obs_names = anndata.obs_names
        obs = pd.DataFrame.from_dict({
            "section": anndata.obs[dataset_section_column],
        })
        spatial = anndata.obsm["spatial"]

        standardized_anndata = ad.AnnData(
            raw_counts,
            obs=obs,
            obsm={"spatial": spatial},
        )
        standardized_anndata.obs["dataset"] = dataset_name
        standardized_anndata.var_names = var_names
        standardized_anndata.obs_names = obs_names

        section = anndata.obs[dataset_section_column].tolist()[0]
        return dataset_name, section, standardized_anndata

    datasets = []

    slide_files = os.listdir("data/seqfish_mouse_organogenesis_imputed/")
    for i, slide in enumerate(slide_files):
        print(f"Processing slide {i + 1} of {len(slide_files)}: {slide}")
        anndata = ad.read_h5ad(f"data/seqfish_mouse_organogenesis_imputed/{slide}")
        dataset_name, section, standardized_anndata = standardize_anndata("seqfish", ["SETD1A", "Setd1a"], anndata)
        datasets.append((dataset_name, section, standardized_anndata))

    return datasets


datasets = load()

os.makedirs("mouse-organogenesis-transcriptomics", exist_ok=True)
for dataset in datasets:
    dataset_name, section, standardized_anndata = dataset
    file_name = f"{dataset_name}__{section}.h5ad"

    local_file = os.path.join(
        "mouse-organogenesis-transcriptomics", file_name
    )
    standardized_anndata.write(local_file)
