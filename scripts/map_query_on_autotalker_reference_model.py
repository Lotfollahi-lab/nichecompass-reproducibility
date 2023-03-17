#!/usr/bin/env python
# coding: utf-8

###############################################################################
# Autotalker Query Mapping on Reference Model #
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
from datetime import datetime

import anndata as ad
import mlflow
import scipy.sparse as sp
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
    help="Input dataset name. The adata file name has to be f'{dataset}.h5ad' "
         "if `query_batches` is `None`. If `query_batches` is not `None`, the "
         "adata file names have to be f'{dataset}_{batch}.h5ad' for each batch "
         "in `query_batches`.")
parser.add_argument(
    "--query_batches",
    nargs='+',
    default=None,
    help="Batches of the input dataset used as query. If not `None`, the adata "
         "file names have to be f'{dataset}_{batch}.h5ad' for each batch in "
         "`query_batches`.")
parser.add_argument(
    "--n_neighbors",
    type=int,
    default=12,
    help="Number of neighbors used to compute the spatial neighborhood graphs.")
parser.add_argument(
    "--spatial_key",
    type=str,
    default="spatial",
    help="Key in `adata.obsm` where the spatial coordinates are stored.")
parser.add_argument(
    "--mapping_entity_key",
    type=str,
    default="mapping_entity",
    help="Key in `adata.obsm` where the mapping entities are stored.")

# Gene program mask
parser.add_argument(
    "--nichenet_keep_target_genes_ratio",
    type=float,
    default=0.01,
    help="Ratio how many of the overall top scored target genes are kept.")
parser.add_argument(
    "--nichenet_max_n_target_genes_per_gp",
    type=int,
    default=50,
    help="After this number of genes the genes are clipped from the gp.")
parser.add_argument(
    "--include_mebocost_gps",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether to include mebocost gene programs.")
parser.add_argument(
    "--mebocost_species",
    type=str,
    default="human",
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

# Model
parser.add_argument(
    "--reference_model_label",
    type=str,
    default="reference",
    help="Label of the reference model under which it was saved.")
parser.add_argument(
    "--load_timestamp",
    type=str,
    help="Timestamp of the reference model training run.")
parser.add_argument(
    "--model_label",
    type=str,
    default="query",
    help="Label of the model under which it will be saved.")
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
    nargs='+',
    default=["gene_expr_decoder",
             "graph_decoder"],
    help="s. Autotalker class signature")
parser.add_argument(
    "--n_cond_embed",
    type=int,
    default=180,
    help="s. Autotalker class signature")
parser.add_argument(
    "--log_variational",
    action=argparse.BooleanOptionalAction,
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
    default=0.01,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--lambda_l1_masked",
    type=float,
    default=0.01,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--edge_batch_size",
    type=int,
    default=128,
    help="s. Autotalker train method signature")
parser.add_argument(
    "--node_batch_size",
    type=int,
    default=16,
    help="s. Autotalker train method signature")

args = parser.parse_args()
if args.query_batches is not None:
    query_batches = [batch for batch in args.query_batches]

# Get time of script execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")

print(f"Run timestamp: {current_timestamp}.")
print("Script arguments:")
print(sys.argv)

# Set mlflow experiment
experiment = mlflow.set_experiment(f"{args.dataset}_{args.model_label}")
mlflow_experiment_id = experiment.experiment_id
mlflow.log_param("timestamp", current_timestamp)

###############################################################################
### 1.3 Configure Paths and Create Directories ###
###############################################################################

model_artifacts_folder_path = f"../artifacts/{args.dataset}/models/" \
                                   f"{args.load_timestamp}"
gp_data_folder_path = "../datasets/gp_data" # gene program data
srt_data_folder_path = "../datasets/srt_data" # spatially-resolved
                                              # transcriptomics data
srt_data_gold_folder_path = f"{srt_data_folder_path}/gold"
nichenet_ligand_target_mx_file_path = gp_data_folder_path + \
                                      "/nichenet_ligand_target_matrix.csv"
omnipath_lr_interactions_file_path = gp_data_folder_path + \
                                     "/omnipath_lr_interactions.csv"
os.makedirs(model_artifacts_folder_path, exist_ok=True)

###############################################################################
## 2. Gene Program Mask ##
###############################################################################

print("\nPreparing the gene program mask...")
# OmniPath gene programs
omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
    min_curation_effort=0,
    load_from_disk=True,
    save_to_disk=False,
    file_path=omnipath_lr_interactions_file_path,
    plot_gp_gene_count_distributions=False)

# NicheNet gene programs
nichenet_gp_dict = extract_gp_dict_from_nichenet_ligand_target_mx(
    keep_target_genes_ratio=args.nichenet_keep_target_genes_ratio,
    max_n_target_genes_per_gp=args.nichenet_max_n_target_genes_per_gp,
    load_from_disk=True,
    save_to_disk=False,
    file_path=nichenet_ligand_target_mx_file_path,
    plot_gp_gene_count_distributions=False)

# Combine gene programs into one dictionary
combined_gp_dict = dict(omnipath_gp_dict)
combined_gp_dict.update(nichenet_gp_dict)

# Mebocost gene programs
if args.include_mebocost_gps:
    mebocost_gp_dict = extract_gp_dict_from_mebocost_es_interactions(
    dir_path=f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps/",
    species=args.mebocost_species,
    genes_uppercase=True,
    plot_gp_gene_count_distributions=False)

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
## 3. Data ##
###############################################################################

###############################################################################
### 3.1 Load Data & Compute Spatial Neighbor Graph ###
###############################################################################

adata_batch_list = []
for batch in query_batches:
    print(f"Processing batch {batch}...")
    print(f"\nLoading data...")
    adata = ad.read_h5ad(
        f"{srt_data_gold_folder_path}/{args.dataset}_{batch}.h5ad")
    adata.obs[args.mapping_entity_key] = "query"

    print(f"\nComputing spatial neighborhood graph...")
    # Compute (separate) spatial neighborhood graphs
    sq.gr.spatial_neighbors(adata,
                            coord_type="generic",
                            spatial_key=args.spatial_key,
                            n_neighs=args.n_neighbors)
    # Make adjacency matrix symmetric
    adata.obsp[args.adj_key] = (
        adata.obsp[args.adj_key].maximum(
            adata.obsp[args.adj_key].T))

    adata_batch_list.append(adata)

adata_query = ad.concat(adata_batch_list, join="inner")

# Combine spatial neighborhood graphs as disconnected components
batch_connectivities = []
len_before_batch = 0
for i in range(len(adata_batch_list)):
    if i == 0: # first batch
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[0].shape[0],
             (adata_query.shape[0] -
              adata_batch_list[0].shape[0])))
        batch_connectivities.append(sp.hstack(
            (adata_batch_list[0].obsp[args.adj_key],
             after_batch_connectivities_extension)))
    elif i == (len(adata_batch_list) - 1): # last batch
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
             (adata_query.shape[0] -
              adata_batch_list[i].shape[0])))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
             adata_batch_list[i].obsp[args.adj_key])))
    else: # middle batches
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0], len_before_batch))
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
             (adata_query.shape[0] -
              adata_batch_list[i].shape[0] -
              len_before_batch)))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
             adata_batch_list[i].obsp[args.adj_key],
             after_batch_connectivities_extension)))
    len_before_batch += adata_batch_list[i].shape[0]
        
connectivities = sp.vstack(batch_connectivities)
adata_query.obsp[args.adj_key] = connectivities

###############################################################################
### 3.2 Add Gene Program Mask to Data ###
###############################################################################

# Add the gene program dictionary as binary masks to the adata for model 
# training
add_gps_from_gp_dict_to_adata(
    gp_dict=combined_new_gp_dict,
    adata=adata_query,
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
n_hidden_encoder = len(adata_query.uns[args.gp_names_key])

###############################################################################
## 4. Initialize, Train & Save Model ##
###############################################################################

print("\nTraining model...")
# Load model trained on reference data for transfer learning with query data
# Freeze all weights except for conditional weights
model = Autotalker.load(
    dir_path=model_artifacts_folder_path + f"/{args.reference_model_label}",
    adata=adata_query,
    adata_file_name=f"{args.dataset}_reference.h5ad",
    gp_names_key=args.gp_names_key,
    unfreeze_all_weights=False,
    unfreeze_cond_embed_weights=True)

# Train model
model.train(n_epochs=args.n_epochs,
            n_epochs_all_gps=args.n_epochs_all_gps,
            lr=args.lr,
            lambda_edge_recon=args.lambda_edge_recon,
            lambda_gene_expr_recon=args.lambda_gene_expr_recon,
            lambda_group_lasso=args.lambda_group_lasso,
            lambda_l1_masked=args.lambda_l1_masked,
            edge_batch_size=args.edge_batch_size,
            node_batch_size=args.node_batch_size,
            mlflow_experiment_id=mlflow_experiment_id,
            verbose=False)

print("\nSaving model...")
# Save trained model
model.save(dir_path=model_artifacts_folder_path + f"/{args.model_label}",
           overwrite=True,
           save_adata=True,
           adata_file_name=f"{args.dataset}_query.h5ad")
