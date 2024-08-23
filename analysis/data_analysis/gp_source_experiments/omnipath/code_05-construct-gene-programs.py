import pickle
import os

from nichecompass.utils import (extract_gp_dict_from_mebocost_es_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                filter_and_combine_gp_dict_gps)


config = {
    "sources": ["omnipath"],
    "omnipath": {
        "min_curation_effort": 0,
    },
    "nichenet": {"keep_target_genes_ratio": 1.,
                 "max_n_target_genes_per_gp": 250,
    },
    "combine_gps": {
        "overlap_thresh_source_genes": 0.9,
        "overlap_thresh_target_genes": 0.9,
        "overlap_thresh_genes": 0.9,
    }
}


def generate_gene_programs(gene_orthologs_mapping_file_path, mebocost_enzyme_sensor_interactions_folder_path, config):

    # Retrieve OmniPath GPs (source: ligand genes; target: receptor genes)
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species="mouse",
        min_curation_effort=config["omnipath"]["min_curation_effort"],
        load_from_disk=False,
        save_to_disk=False,
        gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
        plot_gp_gene_count_distributions=False
    )

    # Retrieve MEBOCOST GPs (source: enzyme genes; target: sensor genes)
    mebocost_gp_dict = extract_gp_dict_from_mebocost_es_interactions(
        dir_path=mebocost_enzyme_sensor_interactions_folder_path,
        species="mouse",
        plot_gp_gene_count_distributions=False
    )

    # Retrieve NicheNet GPs (source: ligand genes; target: receptor genes, target genes)
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species="mouse",
        version="v2",
        keep_target_genes_ratio=config["nichenet"]["keep_target_genes_ratio"],
        max_n_target_genes_per_gp=config["nichenet"]["max_n_target_genes_per_gp"],
        load_from_disk=False,
        save_to_disk=False,
        plot_gp_gene_count_distributions=False
    )

    # Add GPs into one combined dictionary for model training
    combined_gp_dict = {}
    if "omnipath" in config["sources"]:
        combined_gp_dict.update(omnipath_gp_dict)
    if "mebocost" in config["sources"]:
        combined_gp_dict.update(mebocost_gp_dict)
    if "nichenet" in config["sources"]:
        combined_gp_dict.update(nichenet_gp_dict)

    # Filter and combine GPs to avoid overlaps
    combined_new_gp_dict = filter_and_combine_gp_dict_gps(
        gp_dict=combined_gp_dict,
        gp_filter_mode="subset",
        combine_overlap_gps=True,
        overlap_thresh_source_genes=config["combine_gps"]["overlap_thresh_source_genes"],
        overlap_thresh_target_genes=config["combine_gps"]["overlap_thresh_target_genes"],
        overlap_thresh_genes=config["combine_gps"]["overlap_thresh_genes"],)

    print("Number of gene programs before filtering and combining: "
          f"{len(combined_gp_dict)}.")
    print(f"Number of gene programs after filtering and combining: "
          f"{len(combined_new_gp_dict)}.")

    return combined_new_gp_dict



mebocost_dir = "mebocost-gene-programs"

mouse_gene_orthologs_dir = "mouse-gene-orthologs"

gene_programs = generate_gene_programs(
    config=config,
    gene_orthologs_mapping_file_path=os.path.join(mouse_gene_orthologs_dir, "human_mouse_gene_orthologs.csv"),
    mebocost_enzyme_sensor_interactions_folder_path=mebocost_dir
)

os.makedirs("processed-gene-programs:omnipath", exist_ok=True)

file_name = "gene_programs.pkl"

local_file = os.path.join("processed-gene-programs:omnipath", file_name)
with open(local_file, "wb") as file:
    pickle.dump(gene_programs, file)
