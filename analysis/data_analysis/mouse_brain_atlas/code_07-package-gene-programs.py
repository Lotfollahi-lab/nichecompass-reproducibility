import os
import pathlib
import anndata as ad
import pickle

from nichecompass.utils import add_gps_from_gp_dict_to_adata


def package_gene_programs(adata, gene_program_dict):
    gp_targets_mask_key = "nichecompass_gp_targets"
    gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
    gp_sources_mask_key = "nichecompass_gp_sources"
    gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
    gp_names_key = "nichecompass_gp_names"

    packaged_anndata = adata.copy()

    # Add the GP dictionary as binary masks to the adata
    add_gps_from_gp_dict_to_adata(
        gp_dict=gene_program_dict,
        adata=packaged_anndata,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        gp_names_key=gp_names_key,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None)
    return packaged_anndata


dataset = "dataset-with-neighborhood-graph"
anndata = ad.read_h5ad(os.path.join(dataset, "dataset_with_neighborhood_graph.h5ad"))
gene_programs_directory = "processed-gene-programs"
with open(os.path.join(gene_programs_directory, "gene_programs.pkl"), "rb") as file:
    gene_program_dict = pickle.load(file)
packaged_anndata = package_gene_programs(anndata, gene_program_dict)
os.makedirs("dataset-with-gene-programs", exist_ok=True)
file_name = "dataset_with_gene_programs.h5ad"
local_file = os.path.join("dataset-with-gene-programs", file_name)
packaged_anndata.write_h5ad(pathlib.Path(local_file))
