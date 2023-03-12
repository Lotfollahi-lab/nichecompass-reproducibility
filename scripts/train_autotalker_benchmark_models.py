#!/usr/bin/env python
# coding: utf-8

###############################################################################
# Autotalker Benchmark Models Training #
###############################################################################

###############################################################################
## 1. Setup ##
###############################################################################

###############################################################################
### 1.1 Import Libraries ###
###############################################################################

import sys
sys.path.append("../../autotalker")

import argparse
import os
import sys
import time
from datetime import datetime

import anndata as ad
import mlflow
import squidpy as sq

from autotalker.models import Autotalker
from autotalker.utils import (add_gps_from_gp_dict_to_adata,
                              extract_gp_dict_from_mebocost_es_interactions,
                              extract_gp_dict_from_nichenet_ligand_target_mx,
                              extract_gp_dict_from_omnipath_lr_interactions,
                              filter_and_combine_gp_dict_gps)

###############################################################################
### 1.2 Define Parameters ###
###############################################################################

parser = argparse.ArgumentParser(description=os.path.basename(__file__))

# Dataset
parser.add_argument(
    "--dataset",
    type=str,
    help="Input dataset. The adata file name has to be f'{dataset}.h5ad' ")
parser.add_argument(
    "--adata_new_name",
    type=str,
    default=None,
    help="Name of adata to which to append results. If `None`, create.")
parser.add_argument(
    "--n_neighbors_list",
    nargs="+",
    default="4 4 8 8 12 12 16 16 20 20",
    help="Number of neighbors per model training run.")
parser.add_argument(
    "--cell_type_key",
    type=str,
    default="cell_type",
    help="Key in `adata.obs` where the cell types are stored.")
parser.add_argument(
    "--spatial_key",
    type=str,
    default="spatial",
    help="Key in `adata.obsm` where the spatial coordinates are stored.")

# Gene program mask
parser.add_argument(
    "--nichenet_keep_target_genes_ratio",
    type=float,
    default=0.01,
    help="Ratio how many of the overall top scored target genes are kept.")
parser.add_argument(
    "--nichenet_max_n_target_genes_per_gp",
    type=int,
    default=20000,
    help="After this number of genes the genes are clipped from the gp.")
parser.add_argument(
    "--include_mebocost_gps",
    type=bool,
    default=True,
    help="Indicator whether to include mebocost gene programs.")
parser.add_argument(
    "--mebocost_species",
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
    type=bool,
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

# Model
parser.add_argument(
    "--model_name",
    type=str,
    default="autotalker",
    help="Name of the model.")
parser.add_argument(
    "--edge_batch_size_list",
    nargs="+",
    default="512 512 256 256 256 256 128 128 128 128",
    help="Edge batch sizes per model training run.")
parser.add_argument(
    "--node_batch_size_list",
    nargs="+",
    default="64 64 32 32 32 32 16 16 16 16",
    help="Node batch sizes per model training run.")
parser.add_argument(
    "--counts_key",
    type=str,
    default="counts",
    help="s. Autotalker class signature")
parser.add_argument(
    "--adj_key",
    type=str,
    default="spatial_connectivities",
    help="s. Autotalker class signature")
parser.add_argument(
    "--gp_names_key",
    type=str,
    default="autotalker_gp_names",
    help="s. Autotalker class signature")
parser.add_argument(
    "--active_gp_names_key",
    type=str,
    default="autotalker_active_gp_names",
    help="s. Autotalker class signature")
parser.add_argument(
    "--gp_targets_mask_key",
    type=str,
    default="autotalker_gp_targets",
    help="s. Autotalker class signature")
parser.add_argument(
    "--gp_sources_mask_key",
    type=str,
    default="autotalker_gp_sources",
    help="s. Autotalker class signature")
parser.add_argument(
    "--latent_key",
    type=str,
    default="autotalker_latent",
    help="s. Autotalker class signature")
parser.add_argument(
    "--condition_key",
    type=str,
    default="batch",
    help="s. Autotalker class signature")
parser.add_argument(
    "--active_gp_thresh_ratio",
    type=float,
    default=0.03,
    help="s. Autotalker class signature")
parser.add_argument(
    "--gene_expr_recon_dist",
    type=str,
    default="nb",
    help="s. Autotalker class signature")
parser.add_argument(
    "--cond_embed_injection",
    type=list,
    default=["gene_expr_decoder",
             "graph_decoder"],
    help="s. Autotalker class signature")
parser.add_argument(
    "--n_cond_embed",
    type=int,
    default=20,
    help="s. Autotalker class signature")
parser.add_argument(
    "--log_variational",
    type=bool,
    default=True,
    help="s. Autotalker class signature") # counts as input
parser.add_argument(
    "--n_layers_encoder",
    type=int,
    default=1,
    help="s. Autotalker class signature")
parser.add_argument(
    "--conv_layer_encoder",
    type=str,
    default="gcnconv",
    help="s. Autotalker class signature")
parser.add_argument(
    "--n_epochs",
    type=int,
    default=40,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--n_epochs_all_gps",
    type=int,
    default=20,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--lr",
    type=float,
    default=0.001,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--lambda_edge_recon",
    type=float,
    default=1000.,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--lambda_gene_expr_recon",
    type=float,
    default=1.,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--lambda_group_lasso",
    type=float,
    default=0.,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--lambda_l1_masked",
    type=float,
    default=0.,
    help="s. Autotalker train method signature")

# Other
parser.add_argument(
    "--seeds",
    nargs="+",
    default="0 1 2 3 4 5 6 7 8 9",
    help="Random seeds per model training run.")
parser.add_argument(
    "--run_index",
    nargs="+",
    default="1 2 3 4 5 6 7 8 9 10",
    help="Index per model training run.")

args = parser.parse_args()
n_neighbors_list = [int(n_neighbors) for n_neighbors in args.n_neighbors_list]
edge_batch_size_list = [
    int(edge_batch_size) for edge_batch_size in args.edge_batch_size_list]
node_batch_size_list = [
    int(node_batch_size) for node_batch_size in args.node_batch_size_list]
seeds = [int(seed) for seed in args.seeds]
run_index = [int(run_idx) for run_idx in args.run_index]

print(n_neighbors_list)

# Get time of script execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")

print(f"Run timestamp: {current_timestamp}.")
print("Script arguments:")
print(sys.argv)

###############################################################################
### 1.3 Configure Paths and Create Directories ###
###############################################################################

figure_folder_path = f"../../figures/{args.dataset}/method_benchmarking/" \
                             f"{args.model_name}/{current_timestamp}"
gp_data_folder_path = "../datasets/gp_data" # gene program data
srt_data_folder_path = "../datasets/srt_data" # spatially-resolved
                                              # transcriptomics data
srt_data_gold_folder_path = f"{srt_data_folder_path}/gold"
nichenet_ligand_target_mx_file_path = gp_data_folder_path + \
                                      "/nichenet_ligand_target_matrix.csv"
omnipath_lr_interactions_file_path = gp_data_folder_path + \
                                     "/omnipath_lr_interactions.csv"
os.makedirs(figure_folder_path, exist_ok=True)

###############################################################################
## 2. Prepare Gene Program Mask ##
###############################################################################

print("\nPreparing the gene program mask...")
# NicheNet gene programs
nichenet_gp_dict = extract_gp_dict_from_nichenet_ligand_target_mx(
    keep_target_genes_ratio=args.nichenet_keep_target_genes_ratio,
    max_n_target_genes_per_gp=args.nichenet_max_n_target_genes_per_gp,
    load_from_disk=True,
    save_to_disk=False,
    file_path=nichenet_ligand_target_mx_file_path,
    plot_gp_gene_count_distributions=False)

# OmniPath gene programs
omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
    min_curation_effort=0,
    load_from_disk=True,
    save_to_disk=False,
    file_path=omnipath_lr_interactions_file_path,
    plot_gp_gene_count_distributions=False)

# Mebocost gene programs
if args.include_mebocost_gps:
    mebocost_gp_dict = extract_gp_dict_from_mebocost_es_interactions(
    dir_path=f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps/",
    species=args.mebocost_species,
    genes_uppercase=True,
    plot_gp_gene_count_distributions=False)
    
# Combine gene programs into one dictionary
combined_gp_dict = dict(nichenet_gp_dict)
combined_gp_dict.update(omnipath_gp_dict)
if args.include_mebocost_gps:
    combined_gp_dict.update(mebocost_gp_dict)
    
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
## 3. Train Models ##
###############################################################################

if args.adata_new_name is None:
    # Create new adata to store benchmarking results
    adata_original = ad.read_h5ad(
        f"{srt_data_gold_folder_path}/{args.dataset}.h5ad")
    adata_new = ad.AnnData(sp.csr_matrix(
        (adata_original.shape[0], adata_original.shape[1]),
        dtype=np.float32))
    adata_new.obs["cell_type"] = adata_original.obs[args.cell_type_key].values
    del(adata_original)
else:
    adata_new = ad.read_h5ad(
        f"{srt_data_gold_folder_path}/{args.adata_new_name}.h5ad")
    
for run_number, n_neighbors in zip(run_index, n_neighbors_list):
    # Load data
    adata = ad.read_h5ad(f"{srt_data_gold_folder_path}/{args.dataset}.h5ad")
    adata.obs[args.condition_key] == "batch1"

    # Add the gene program dictionary as binary masks to the adata for model
    # training
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

    # Determine dimensionality of hidden encoder (in case n_layers_encoder > 1)
    n_hidden_encoder = len(adata.uns[args.gp_names_key])
    
    # Compute spatial neighborhood graph
    sq.gr.spatial_neighbors(adata,
                            coord_type="generic",
                            spatial_key=args.spatial_key,
                            n_neighs=n_neighbors)

    # Make adjacency matrix symmetric
    adata.obsp["spatial_connectivities"] = (
        adata.obsp["spatial_connectivities"].maximum(
            adata.obsp["spatial_connectivities"].T))

    start_time = time.time()
    
    print("\nTraining model...")
    # Initialize model
    model = Autotalker(adata,
                       counts_key=args.counts_key,
                       adj_key=args.adj_key,
                       condition_key=args.condition_key,
                       cond_embed_injection=args.cond_embed_injection,
                       n_cond_embed=args.n_cond_embed,
                       gp_names_key=args.gp_names_key,
                       active_gp_names_key=args.active_gp_names_key,
                       gp_targets_mask_key=args.gp_targets_mask_key,
                       gp_sources_mask_key=args.gp_sources_mask_key,
                       latent_key=args.latent_key,
                       active_gp_thresh_ratio=args.active_gp_thresh_ratio,
                       gene_expr_recon_dist=args.gene_expr_recon_dist,
                       n_layers_encoder=args.n_layers_encoder,
                       conv_layer_encoder=args.conv_layer_encoder,
                       n_hidden_encoder=n_hidden_encoder,
                       log_variational=args.log_variational)

    # Train model
    model.train(n_epochs=args.n_epochs,
                n_epochs_all_gps=args.n_epochs_all_gps,
                lr=args.lr,
                lambda_edge_recon=args.lambda_edge_recon,
                lambda_gene_expr_recon=args.lambda_gene_expr_recon,
                lambda_group_lasso=args.lambda_group_lasso,
                lambda_l1_masked=args.lambda_l1_masked,
                edge_batch_size=edge_batch_size_list[run_number-1],
                node_batch_size=node_batch_size_list[run_number-1],
                seed=seeds[run_number-1],
                verbose=True)
    
    # Measure time for model training
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"Duration of model training in run {run_number}: {int(hours)} "
          f"hours, {int(minutes)} minutes and {int(seconds)} seconds.")
    
    # Store model training duration
    adata_new.uns[
        f"{args.model_name}_model_training_duration_run{run_number}"] = (
        elapsed_time)

    # Store latent representation
    adata_new.obsm[args.latent_key + f"_run{run_number}"] = (
        adata.obsm[args.latent_key])

    # Store intermediate adata to disk
    adata_new.write(
        f"{srt_data_gold_folder_path}/{args.dataset}_{args.model_name}.h5ad")

# Store final adata to disk
adata_new.write(
    f"{srt_data_gold_folder_path}/{args.dataset}_{args.model_name}.h5ad")