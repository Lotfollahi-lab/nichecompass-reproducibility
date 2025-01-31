#!/usr/bin/env python
# coding: utf-8

###############################################################################
# NicheCompass Benchmarking Models Training #
###############################################################################

###############################################################################
## 1. Setup ##
###############################################################################

###############################################################################
### 1.1 Import Libraries ###
###############################################################################

import argparse
import gc
import os
import sys
import time
from datetime import datetime

import anndata as ad
import mlflow
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import squidpy as sq
import torch

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                extract_gp_dict_from_mebocost_ms_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                filter_and_combine_gp_dict_gps,
                                filter_and_combine_gp_dict_gps_v2,
                                get_unique_genes_from_gp_dict)

###############################################################################
### 1.2 Define Parameters ###
###############################################################################

parser = argparse.ArgumentParser(description=os.path.basename(__file__))

def none_or_value(value):
    if value == "None":
        return None
    return value

def none_or_int(value):
    if value == "None":
        return None
    return int(value)

def none_or_bool(value):
    if value == "None":
        return None
    return value == "True"

# Benchmarking-specific
parser.add_argument(
    "--adata_new_name",
    type=none_or_value,
    default=None,
    help="Name of adata to which to append results. If `None`, create new "
         "one.")
parser.add_argument(
    "--n_neighbors_list",
    nargs="+",
    default="4 4 8 8 12 12 16 16",
    help="Number of neighbors per model training run.")
parser.add_argument(
    "--edge_batch_size_list",
    nargs="+",
    default="256 256 256 256 128 128 64 64",
    help="Edge batch sizes per model training run.")
parser.add_argument(
    "--node_batch_size_list",
    type=none_or_value,
    nargs="+",
    default="None None None None None None None None None None",
    help="Node batch sizes per model training run.")
parser.add_argument(
    "--seeds",
    nargs="+",
    default="0 1 2 3 4 5 6 7",
    help="Random seeds per model training run.")
parser.add_argument(
    "--run_index",
    nargs="+",
    default="1 2 3 4 5 6 7 8",
    help="Index per model training run.")

# Gene program mask
parser.add_argument(
    "--nichenet_keep_target_genes_ratio",
    type=float,
    default=1.,
    help="Ratio how many of the overall top scored target genes are kept.")
parser.add_argument(
    "--nichenet_max_n_target_genes_per_gp",
    type=int,
    default=250,
    help="After this number of genes the genes are clipped from the gp.")
parser.add_argument(
    "--include_mebocost_gps",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="Indicator whether to include mebocost gene programs.")
parser.add_argument(
    "--species",
    type=str,
    default="mouse",
    help="Species that is used for the retrieval of mebocost gene programs.")
parser.add_argument(
    "--gp_filter_mode",
    type=str,
    default="subset",
    help="Which kind of gene programs are filtered.")
parser.add_argument(
    "--combine_overlap_gps",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="Indicator whether to combine overlapping gene programs.")
parser.add_argument(
    "--overlap_thresh_source_genes",
    type=float,
    default=0.9,
    help="Threshold for source genes above which gene programs are combined.")
parser.add_argument(
    "--overlap_thresh_target_genes",
    type=float,
    default=0.9,
    help="Threshold for target genes above which gene programs are combined.")
parser.add_argument(
    "--overlap_thresh_genes",
    type=float,
    default=0.9,
    help="Threshold for overall genes above which gene programs are combined.")

# Data
parser.add_argument(
    "--dataset",
    type=str,
    help="Input dataset name. The adata file name has to be f'{dataset}.h5ad' "
         "if `reference_batches` is `None`. If `reference_batches` is not "
         "`None`, the adata file names have to be f'{dataset}_{batch}.h5ad' "
         "for each batch in `reference_batches`.")
parser.add_argument(
    "--reference_batches",
    nargs='+',
    type=none_or_value,
    default=None,
    help="Batches of the input dataset used as reference. If not `None`, the "
         "adata file names have to be f'{dataset}_{batch}.h5ad' for each batch"
         " in `reference_batches`.")
parser.add_argument(
    "--counts_key",
    type=none_or_value,
    default="counts",
    help="s. NicheCompass class signature.")
parser.add_argument(
    "--condition_key",
    type=str,
    default="batch",
    help="s. NicheCompass class signature.")
parser.add_argument(
    "--cat_covariates_keys",
    nargs='+',
    type=none_or_value,
    default=None,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--cat_covariates_no_edges",
    nargs='+',
    type=none_or_bool,
    default=None,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--spatial_key",
    type=str,
    default="spatial",
    help="Key in `adata.obsm` where the spatial coordinates are stored.")
parser.add_argument(
    "--adj_key",
    type=str,
    default="spatial_connectivities",
    help="s. NicheCompass class signature.")
parser.add_argument(
    "--cell_type_key",
    type=str,
    default="cell_type",
    help="Key in `adata.obs` where the cell types are stored.")
parser.add_argument(
    "--mapping_entity_key",
    type=str,
    default="mapping_entity",
    help="Key in `adata.obsm` where the mapping entities will be stored.")
parser.add_argument(
    "--filter_genes",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether genes should be filtered.")
parser.add_argument(
    "--n_hvg",
    type=int,
    default=0,
    help="Number of highly variable genes that are kept if `filter_genes` is "
         "`True`.")
parser.add_argument(
    "--n_svg",
    type=int,
    default=3000,
    help="Number of spatially variable genes that are kept if `filter_genes` "
         "is `True`.")
parser.add_argument(
    "--gp_targets_mask_key",
    type=str,
    default="nichecompass_gp_targets",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--gp_sources_mask_key",
    type=str,
    default="nichecompass_gp_sources",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--gp_names_key",
    type=str,
    default="nichecompass_gp_names",
    help="s. NicheCompass class signature")

# Model
parser.add_argument(
    "--model_label",
    type=str,
    default="single_sample_method_benchmarking",
    help="Label of the models under which they will be saved.")
parser.add_argument(
    "--active_gp_names_key",
    type=str,
    default="nichecompass_active_gp_names",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--latent_key",
    type=str,
    default="nichecompass_latent",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--n_addon_gp",
    type=int,
    default=10,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--active_gp_thresh_ratio",
    type=float,
    default=0.05,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--gene_expr_recon_dist",
    type=str,
    default="nb",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--cat_covariates_embeds_injection",
    nargs='+',
    default=["gene_expr_decoder"],
    help="s. NicheCompass class signature")
parser.add_argument(
    "--cat_covariates_embeds_nums",
    nargs='+',
    type=none_or_int,
    default=None,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--log_variational",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="s. NicheCompass class signature") # counts as input
parser.add_argument(
    "--node_label_method",
    type=str,
    default="one-hop-norm",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--n_layers_encoder",
    type=int,
    default=1,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--n_fc_layers_encoder",
    type=int,
    default=1,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--conv_layer_encoder",
    type=str,
    default="gcnconv",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--n_hidden_encoder",
    type=none_or_int,
    default=None,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--n_epochs",
    type=int,
    default=400,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--n_epochs_all_gps",
    type=int,
    default=25,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--n_epochs_no_cat_covariates_contrastive",
    type=int,
    default=5,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lr",
    type=float,
    default=0.001,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_edge_recon",
    type=float,
    default=5000000.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_gene_expr_recon",
    type=float,
    default=3000.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_cat_covariates_contrastive",
    type=float,
    default=0.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--contrastive_logits_pos_ratio",
    type=float,
    default=0.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--contrastive_logits_neg_ratio",
    type=float,
    default=0.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_group_lasso",
    type=float,
    default=0.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_l1_masked",
    type=float,
    default=0.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_l1_addon",
    type=float,
    default=0.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--n_sampled_neighbors",
    type=int,
    default=-1,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--use_new_gp_mask",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether new gene program mask should be used.")

# Other
parser.add_argument(
    "--timestamp_suffix",
    type=str,
    default="",
    help="Timestamp suffix for saving artifacts if timestamp overlaps")

args = parser.parse_args()

n_neighbors_list = [int(n_neighbors) for n_neighbors in args.n_neighbors_list]
edge_batch_size_list = [
    int(edge_batch_size) for edge_batch_size in args.edge_batch_size_list]
node_batch_size_list = [
    int(node_batch_size) if node_batch_size is not None else None for 
    node_batch_size in args.node_batch_size_list]
seeds = [int(seed) for seed in args.seeds]
run_index = [int(run_idx) for run_idx in args.run_index]

if args.reference_batches == [None]:
    args.reference_batches = None
if args.cat_covariates_embeds_injection == [None]:
    args.cat_covariates_embeds_injection = []
if args.cat_covariates_keys == [None]:
    args.cat_covariates_keys = None
if args.cat_covariates_no_edges == [None]:
    args.cat_covariates_no_edges = None
if args.cat_covariates_embeds_nums == [None]:
    args.cat_covariates_embeds_nums = None

# Get time of script execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")

print(f"Run timestamp: {current_timestamp + args.timestamp_suffix}.")
print("Script arguments:")
print(sys.argv)

###############################################################################
### 1.3 Configure Paths and Create Directories ###
###############################################################################

root_folder_path = os.path.dirname(
    os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
artifacts_folder_path = f"{root_folder_path}/artifacts"
model_folder_path = f"{artifacts_folder_path}/{args.dataset}/models/" \
                    f"{args.model_label}/{current_timestamp}" \
                     f"{args.timestamp_suffix}"
result_folder_path = f"{artifacts_folder_path}/{args.dataset}/results/" \
                     f"{args.model_label}/{current_timestamp}" \
                     f"{args.timestamp_suffix}"
gp_data_folder_path = f"{root_folder_path}/datasets/gp_data" # gene program 
                                                             # data
ga_data_folder_path = f"{root_folder_path}/datasets/ga_data" # gene annotation
                                                             # data
so_data_folder_path = f"{root_folder_path}/datasets/st_data" # spatial omics
                                                              # data
so_data_gold_folder_path = f"{so_data_folder_path}/gold"
nichenet_lr_network_file_path = gp_data_folder_path + \
                                "/nichenet_lr_network_v2_" \
                                f"{args.species}.csv"
nichenet_ligand_target_matrix_file_path = gp_data_folder_path + \
                                          "/nichenet_ligand_target_matrix_" \
                                          f"v2_{args.species}.csv"
omnipath_lr_network_file_path = gp_data_folder_path + \
                                     "/omnipath_lr_network.csv"
gene_orthologs_mapping_file_path = ga_data_folder_path + \
                                   "/human_mouse_gene_orthologs.csv"
os.makedirs(model_folder_path, exist_ok=True)
os.makedirs(result_folder_path, exist_ok=True)

###############################################################################
## 2. Prepare Gene Program Mask ##
###############################################################################

print("\nPreparing the gene program mask...")
if args.use_new_gp_mask:
    min_curation_effort = 2
else:
    min_curation_effort = 0

# OmniPath gene programs
omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
    species=args.species,
    min_curation_effort=min_curation_effort,
    load_from_disk=True,
    save_to_disk=False,
    lr_network_file_path=omnipath_lr_network_file_path,
    gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
    plot_gp_gene_count_distributions=False)

omnipath_genes = get_unique_genes_from_gp_dict(
    gp_dict=omnipath_gp_dict,
    retrieved_gene_entities=["sources", "targets"])

# NicheNet gene programs
nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
    species=args.species,
    version="v2",
    keep_target_genes_ratio=args.nichenet_keep_target_genes_ratio,
    max_n_target_genes_per_gp=args.nichenet_max_n_target_genes_per_gp,
    load_from_disk=True,
    save_to_disk=False,
    lr_network_file_path=nichenet_lr_network_file_path,
    ligand_target_matrix_file_path=nichenet_ligand_target_matrix_file_path,
    gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
    plot_gp_gene_count_distributions=False)

nichenet_lr_genes = get_unique_genes_from_gp_dict(
    gp_dict=nichenet_gp_dict,
    retrieved_gene_categories=["ligand", "receptor"])

# Combine gene programs into one dictionary
if args.use_new_gp_mask:
    gp_dicts = [omnipath_gp_dict, nichenet_gp_dict]
combined_gp_dict = dict(omnipath_gp_dict)
combined_gp_dict.update(nichenet_gp_dict)

if args.filter_genes:
    # Get gene program relevant genes
    gp_relevant_genes = [gene.upper() for gene in list(set(
        omnipath_genes + nichenet_lr_genes))]

# Mebocost gene programs
if args.include_mebocost_gps:
    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
    dir_path=f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps",
    species=args.species,
    plot_gp_gene_count_distributions=False)
    
    mebocost_genes = get_unique_genes_from_gp_dict(
        gp_dict=mebocost_gp_dict,
        retrieved_gene_entities=["sources", "targets"])
    
    if args.use_new_gp_mask:
        gp_dicts.append(mebocost_gp_dict)
    combined_gp_dict.update(mebocost_gp_dict)
    
    #if args.filter_genes:
        # Update gene program relevant genes
        #gp_relevant_genes = list(set(gp_relevant_genes + mebocost_genes))
        
if args.use_new_gp_mask:
    combined_new_gp_dict = filter_and_combine_gp_dict_gps_v2(
        gp_dicts,
        verbose=True)
else:
    # Filter and combine gene programs
    combined_new_gp_dict = filter_and_combine_gp_dict_gps(
        gp_dict=combined_gp_dict,
        gp_filter_mode=args.gp_filter_mode,
        combine_overlap_gps=args.combine_overlap_gps,
        overlap_thresh_source_genes=args.overlap_thresh_source_genes,
        overlap_thresh_target_genes=args.overlap_thresh_target_genes,
        overlap_thresh_genes=args.overlap_thresh_genes,
        verbose=False)

print("Number of gene programs before filtering and combining: "
      f"{len(combined_gp_dict)}.")
print(f"Number of gene programs after filtering and combining: "
      f"{len(combined_new_gp_dict)}.")

###############################################################################
## 3. Data ##
###############################################################################

###############################################################################
### 3.1 Create New AnnData to Store Results
###############################################################################

# Create new adata to store results from training runs in storage-efficient way
if args.adata_new_name is None:  
    adata_batch_list = []
    if args.reference_batches is not None:
        for batch in args.reference_batches:
            adata_batch = ad.read_h5ad(
                f"{so_data_gold_folder_path}/{args.dataset}_{batch}.h5ad")
            adata_batch.obs[args.mapping_entity_key] = "reference"
            adata_batch_list.append(adata_batch)
        adata_original = ad.concat(adata_batch_list, join="inner")
    else:
        adata_original = ad.read_h5ad(
            f"{so_data_gold_folder_path}/{args.dataset}.h5ad")

    adata_new = sc.AnnData(sp.csr_matrix(
        (adata_original.shape[0], adata_original.shape[1]),
        dtype=np.float32))
    adata_new.var_names = adata_original.var_names
    adata_new.obs_names = adata_original.obs_names
    adata_new.obs["cell_type"] = adata_original.obs[args.cell_type_key].values
    adata_new.obsm[args.spatial_key] = adata_original.obsm[args.spatial_key]
    if args.cat_covariates_keys is not None:
        for cat_covariate_key in args.cat_covariates_keys:
            adata_new.obs[cat_covariate_key] = (
                adata_original.obs[cat_covariate_key])
    if args.model_label == "sample_integration_method_benchmarking":
        adata_new.obs[args.mapping_entity_key] = (
            adata_original.obs[args.mapping_entity_key])
    del(adata_original)
else:
    adata_new = ad.read_h5ad(
        f"{result_folder_path}/{args.adata_new_name}")
    
###############################################################################
## 4. Train Models ##
###############################################################################
    
for k, (run_number, n_neighbors) in enumerate(zip(run_index,
                                                  n_neighbors_list)):
    # Load data
    adata_batch_list = []
    if args.reference_batches is not None:
        for batch in args.reference_batches:
            print(f"Processing batch {batch}...")
            print("Loading data...")
            adata_batch = ad.read_h5ad(
                f"{so_data_gold_folder_path}/{args.dataset}_{batch}.h5ad")
            adata_batch.obs[args.mapping_entity_key] = "reference"
            print("Computing spatial neighborhood graph...\n")
            # Compute (separate) spatial neighborhood graphs
            sq.gr.spatial_neighbors(adata_batch,
                                    coord_type="generic",
                                    spatial_key=args.spatial_key,
                                    n_neighs=n_neighbors)
            # Make adjacency matrix symmetric
            adata_batch.obsp[args.adj_key] = (
                adata_batch.obsp[args.adj_key].maximum(
                    adata_batch.obsp[args.adj_key].T))
            adata_batch_list.append(adata_batch)
        adata = ad.concat(adata_batch_list, join="inner")

        # Combine spatial neighborhood graphs as disconnected components
        batch_connectivities = []
        len_before_batch = 0
        for i in range(len(adata_batch_list)):
            if i == 0: # first batch
                after_batch_connectivities_extension = sp.csr_matrix(
                    (adata_batch_list[0].shape[0],
                    (adata.shape[0] -
                    adata_batch_list[0].shape[0])))
                batch_connectivities.append(sp.hstack(
                    (adata_batch_list[0].obsp[args.adj_key],
                    after_batch_connectivities_extension)))
            elif i == (len(adata_batch_list) - 1): # last batch
                before_batch_connectivities_extension = sp.csr_matrix(
                    (adata_batch_list[i].shape[0],
                    (adata.shape[0] -
                    adata_batch_list[i].shape[0])))
                batch_connectivities.append(sp.hstack(
                    (before_batch_connectivities_extension,
                    adata_batch_list[i].obsp[args.adj_key])))
            else: # middle batches
                before_batch_connectivities_extension = sp.csr_matrix(
                    (adata_batch_list[i].shape[0], len_before_batch))
                after_batch_connectivities_extension = sp.csr_matrix(
                    (adata_batch_list[i].shape[0],
                    (adata.shape[0] -
                    adata_batch_list[i].shape[0] -
                    len_before_batch)))
                batch_connectivities.append(sp.hstack(
                    (before_batch_connectivities_extension,
                    adata_batch_list[i].obsp[args.adj_key],
                    after_batch_connectivities_extension)))
            len_before_batch += adata_batch_list[i].shape[0]
        connectivities = sp.vstack(batch_connectivities)
        adata.obsp[args.adj_key] = connectivities
    else:
        adata = ad.read_h5ad(
            f"{so_data_gold_folder_path}/{args.dataset}.h5ad")
        # Compute (separate) spatial neighborhood graphs
        sq.gr.spatial_neighbors(adata,
                                coord_type="generic",
                                spatial_key=args.spatial_key,
                                n_neighs=n_neighbors)
        # Make adjacency matrix symmetric
        adata.obsp[args.adj_key] = (
            adata.obsp[args.adj_key].maximum(
                adata.obsp[args.adj_key].T))

    # Filter genes if specified
    if args.filter_genes:
        print("\nFiltering genes...")
        # Filter genes and only keep gp relevant, the 'n_svg' spatially variable,
        # and 'n_hvg' highly variable genes
        gp_dict_genes = get_unique_genes_from_gp_dict(
            gp_dict=combined_new_gp_dict,
            retrieved_gene_entities=["sources", "targets"])
        print(f"Starting with {len(adata.var_names)} genes.")
        sc.pp.filter_genes(adata,
                           min_cells=0)
        print(f"Keeping {len(adata.var_names)} genes after filtering genes with "
              "expression in 0 cells.")

        if args.counts_key is not None:
            hvg_layer = args.counts_key
            if (adata.layers[args.counts_key].astype(int).astype(np.float32).sum() == 
            adata.layers[args.counts_key].sum()): # raw counts
                hvg_flavor = "seurat_v3"
            else: # log normalized counts
                hvg_flavor = "seurat"
        else:
            hvg_layer = None
            if adata.X.astype(int).astype(np.float32).sum() == adata.X.sum():
            # raw counts
                hvg_flavor = "seurat_v3"
            else: # log normalized counts
                hvg_flavor = "seurat"

        if args.n_hvg != 0:
            sc.pp.highly_variable_genes(
                adata,
                layer=hvg_layer,
                n_top_genes=args.n_hvg,
                flavor=hvg_flavor,
                batch_key=args.condition_key,
                subset=False)
        else:
            adata.var["highly_variable"] = False

        # Filter spatially variable genes
        if args.n_svg != 0:
            sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names)
            sv_genes = adata.uns["moranI"].index[:args.n_svg].tolist()
            adata.var["spatially_variable"] = adata.var_names.isin(sv_genes)
        else:
            adata.var["spatially_variable"] = False
        
        gp_relevant_genes = [] # remove gp relevant genes
        adata.var["gp_relevant"] = (
            adata.var.index.str.upper().isin(gp_relevant_genes))
        adata.var["keep_gene"] = (adata.var["gp_relevant"] | 
                                  adata.var["highly_variable"] |
                                  adata.var["spatially_variable"])
        adata = adata[:, adata.var["keep_gene"] == True]
        print(f"Keeping {len(adata.var_names)} spatially variable, highly "
              "variable or gene program relevant genes.")
        #adata = (adata[:, adata.var_names[adata.var_names.str.upper().isin(
    #            [gene.upper() for gene in gp_dict_genes])].sort_values()])
    #print(f"Keeping {len(adata.var_names)} genes after filtering genes not in "
    #      "gp dict.")

    # Add the gene program dictionary as binary masks to the adata for
    # model training
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_new_gp_dict,
        adata=adata,
        genes_uppercase=True,
        gp_targets_mask_key=args.gp_targets_mask_key,
        gp_sources_mask_key=args.gp_sources_mask_key,
        gp_names_key=args.gp_names_key,
        min_genes_per_gp=1,
        min_source_genes_per_gp=0,
        min_target_genes_per_gp=0,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None,
        filter_genes_not_in_masks=False)
    
    # Set mlflow experiment
    experiment = mlflow.set_experiment(f"{args.dataset}_{args.model_label}")
    mlflow_experiment_id = experiment.experiment_id
    mlflow.log_param("timestamp", current_timestamp + args.timestamp_suffix)

    # Log dataset params
    mlflow.log_param("dataset", args.dataset)
    mlflow.log_param("run_number", run_number)
    mlflow.log_param("n_neighbors", n_neighbors)
    
    # Log gp mask params
    mlflow.log_param("nichenet_keep_target_genes_ratio",
                     args.nichenet_keep_target_genes_ratio)
    mlflow.log_param("nichenet_max_n_target_genes_per_gp",
                     args.nichenet_max_n_target_genes_per_gp)
    mlflow.log_param("include_mebocost_gps", args.include_mebocost_gps)
    mlflow.log_param("species", args.species)
    mlflow.log_param("gp_filter_mode", args.gp_filter_mode)
    mlflow.log_param("combine_overlap_gps", args.combine_overlap_gps)
    mlflow.log_param("overlap_thresh_source_genes",
                     args.overlap_thresh_source_genes)
    mlflow.log_param("overlap_thresh_target_genes",
                     args.overlap_thresh_target_genes)
    mlflow.log_param("overlap_thresh_genes", args.overlap_thresh_genes)

    start_time = time.time()

    print("\nTraining model...")
    # Initialize model
    model = NicheCompass(adata,
                         counts_key=args.counts_key,
                         adj_key=args.adj_key,
                         cat_covariates_embeds_injection=args.cat_covariates_embeds_injection,
                         cat_covariates_keys=args.cat_covariates_keys,
                         cat_covariates_no_edges=args.cat_covariates_no_edges,
                         cat_covariates_embeds_nums=args.cat_covariates_embeds_nums,
                         gp_names_key=args.gp_names_key,
                         active_gp_names_key=args.active_gp_names_key,
                         gp_targets_mask_key=args.gp_targets_mask_key,
                         gp_sources_mask_key=args.gp_sources_mask_key,
                         latent_key=args.latent_key,
                         n_addon_gp=args.n_addon_gp,
                         active_gp_thresh_ratio=args.active_gp_thresh_ratio,
                         gene_expr_recon_dist=args.gene_expr_recon_dist,
                         n_fc_layers_encoder=args.n_fc_layers_encoder,
                         conv_layer_encoder=args.conv_layer_encoder,
                         n_hidden_encoder=args.n_hidden_encoder,
                         log_variational=args.log_variational,
                         node_label_method=args.node_label_method)

    # Train model
    model.train(n_epochs=args.n_epochs,
                n_epochs_all_gps=args.n_epochs_all_gps,
                n_epochs_no_cat_covariates_contrastive=args.n_epochs_no_cat_covariates_contrastive,
                lr=args.lr,
                lambda_edge_recon=args.lambda_edge_recon,
                lambda_gene_expr_recon=args.lambda_gene_expr_recon,
                lambda_cat_covariates_contrastive=args.lambda_cat_covariates_contrastive,
                contrastive_logits_pos_ratio=args.contrastive_logits_pos_ratio,
                contrastive_logits_neg_ratio=args.contrastive_logits_neg_ratio,
                lambda_group_lasso=args.lambda_group_lasso,
                lambda_l1_masked=args.lambda_l1_masked,
                lambda_l1_addon=args.lambda_l1_addon,
                edge_batch_size=edge_batch_size_list[k],
                node_batch_size=node_batch_size_list[k],
                n_sampled_neighbors=args.n_sampled_neighbors,
                mlflow_experiment_id=mlflow_experiment_id,
                seed=seeds[k],
                verbose=True)

    # Measure time for model training
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"Duration of model training in run {run_number}: {int(hours)} "
          f"hours, {int(minutes)} minutes and {int(seconds)} seconds.")
    
    print("\nComputing neighbor graph...")
    # Use latent representation for UMAP generation
    sc.pp.neighbors(model.adata,
                    use_rep=f"{args.latent_key}",
                    key_added=f"{args.latent_key}")

    print("\nComputing UMAP embedding...")
    sc.tl.umap(model.adata,
               neighbors_key=f"{args.latent_key}")
    
    # Store latent representation
    adata_new.obsm[args.latent_key + f"_run{run_number}"] = (
        model.adata.obsm[args.latent_key])
    
    # Store latent nearest neighbor graph
    adata_new.obsp[f"{args.latent_key}_run{run_number}_connectivities"] = (
        model.adata.obsp[f"{args.latent_key}_connectivities"])
    adata_new.obsp[f"{args.latent_key}_run{run_number}_distances"] = (
        model.adata.obsp[f"{args.latent_key}_distances"])

    # Store UMAP features
    adata_new.obsm[f"{args.latent_key}_run{run_number}_X_umap"] = (
        model.adata.obsm["X_umap"])
    adata_new.uns[f"{args.latent_key}_run{run_number}_umap"] = (
        model.adata.uns["umap"])

    # Store model training duration
    adata_new.uns[
        f"nichecompass_model_training_duration_run{run_number}"] = (
        elapsed_time)

    # Store intermediate adata to disk
    adata_new.write(
        f"{result_folder_path}/{args.dataset}_{args.model_label}.h5ad") 

    print("\nSaving model...")
    # Save trained model
    model.save(
        dir_path=f"{model_folder_path}/run{run_number}/",
        overwrite=True,
        save_adata=True,
        adata_file_name=f"{args.dataset}_{args.model_label}.h5ad")
    
    # Free memory
    del(adata)
    del(model)
    gc.collect()
    torch.cuda.empty_cache()
    
    mlflow.end_run()

# Store final adata to disk
adata_new.write(
    f"{result_folder_path}/{args.dataset}_{args.model_label}.h5ad") 