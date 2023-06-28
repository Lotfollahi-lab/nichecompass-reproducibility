#!/usr/bin/env python
# coding: utf-8

###############################################################################
# NicheCompass Query Mapping on Reference Model #
###############################################################################

###############################################################################
## 1. Setup ##
###############################################################################

###############################################################################
### 1.1 Import Libraries ###
###############################################################################

import argparse
import os
import sys

import anndata as ad
import mlflow
import scanpy as sc
import scipy.sparse as sp
import squidpy as sq

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                extract_gp_dict_from_mebocost_es_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                filter_and_combine_gp_dict_gps,
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

# Data
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
    type=none_or_value,
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
    help="Key in `adata.obsm` where the mapping entities will be stored.")
parser.add_argument(
    "--gp_names_key",
    type=str,
    default="nichecompass_gp_names",
    help="s. NicheCompass class signature")

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
    "--query_model_label",
    type=str,
    default="query",
    help="Label of the query model under which it will be saved.")
parser.add_argument(
    "--reference_query_model_label",
    type=str,
    default="reference_query",
    help="Label of the reference query model under which it will be saved.")
parser.add_argument(
    "--n_epochs",
    type=int,
    default=200,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--n_epochs_all_gps",
    type=int,
    default=25,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--n_epochs_no_cond_contrastive",
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
    default=500000.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_gene_expr_recon",
    type=float,
    default=300.,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--lambda_cond_contrastive",
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
    "--edge_batch_size",
    type=int,
    default=256,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--node_batch_size",
    type=none_or_int,
    default=None,
    help="s. NicheCompass train method signature")

args = parser.parse_args()

if args.query_batches == [None]:
    args.query_batches = None

print("Script arguments:")
print(sys.argv)

# Set mlflow experiment
experiment = mlflow.set_experiment(
    f"{args.dataset}_{args.reference_query_model_label}")
mlflow_experiment_id = experiment.experiment_id
mlflow.log_param("timestamp", args.load_timestamp)

###############################################################################
### 1.3 Configure Paths and Create Directories ###
###############################################################################

root_folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
artifacts_folder_path = f"{root_folder_path}/artifacts"
model_folder_path = f"{artifacts_folder_path}/{args.dataset}/models"
gp_data_folder_path = f"{root_folder_path}/datasets/gp_data" # gene program
                                                             # data
so_data_folder_path = f"{root_folder_path}/datasets/srt_data" # spatial omics
                                                               # data
so_data_gold_folder_path = f"{so_data_folder_path}/gold"
so_data_results_folder_path = f"{so_data_folder_path}/results"
nichenet_ligand_target_mx_file_path = gp_data_folder_path + \
                                      "/nichenet_ligand_target_matrix.csv"
omnipath_lr_interactions_file_path = gp_data_folder_path + \
                                     "/omnipath_lr_interactions.csv"
os.makedirs(f"{so_data_results_folder_path}/{args.reference_query_model_label}",
            exist_ok=True)
os.makedirs(f"{so_data_results_folder_path}/{args.query_model_label}",
            exist_ok=True)

###############################################################################
## 2. Reference Model ##
###############################################################################

# Retrieve reference model attributes and adata
print("Retrieving reference model...")
reference_model = NicheCompass.load(
    dir_path=model_folder_path + f"/{args.reference_model_label}/" \
             f"{args.load_timestamp}",
    adata_file_name=f"{args.dataset}_{args.reference_model_label}.h5ad",
    gp_names_key=args.gp_names_key)

reference_adata = reference_model.adata
counts_key = reference_model.counts_key_
condition_key = reference_model.condition_key_
adj_key = reference_model.adj_key_
gp_targets_mask_key = reference_model.gp_targets_mask_key_
gp_sources_mask_key = reference_model.gp_sources_mask_key_
target_genes_idx_key = reference_model.target_genes_idx_key_
source_genes_idx_key = reference_model.source_genes_idx_key_
genes_idx_key = reference_model.genes_idx_key_
latent_key = reference_model.latent_key_

genes = reference_model.adata.var_names

###############################################################################
## 3. Data ##
###############################################################################

###############################################################################
### 3.1 Load Data & Compute Spatial Neighbor Graph ###
###############################################################################

adata_batch_list = []
if args.query_batches is not None:
    for batch in args.query_batches:
        print(f"\nProcessing batch {batch}...")
        print("Loading data...")
        adata_batch = ad.read_h5ad(
            f"{so_data_gold_folder_path}/{args.dataset}_{batch}.h5ad")
        adata_batch.obs[args.mapping_entity_key] = "query"
        print("Computing spatial neighborhood graph...")
        # Compute (separate) spatial neighborhood graphs
        sq.gr.spatial_neighbors(adata_batch,
                                coord_type="generic",
                                spatial_key=args.spatial_key,
                                n_neighs=args.n_neighbors)
        # Make adjacency matrix symmetric
        adata_batch.obsp[adj_key] = (
            adata_batch.obsp[adj_key].maximum(
                adata_batch.obsp[adj_key].T))
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
                (adata_batch_list[0].obsp[adj_key],
                after_batch_connectivities_extension)))
        elif i == (len(adata_batch_list) - 1): # last batch
            before_batch_connectivities_extension = sp.csr_matrix(
                (adata_batch_list[i].shape[0],
                (adata.shape[0] -
                adata_batch_list[i].shape[0])))
            batch_connectivities.append(sp.hstack(
                (before_batch_connectivities_extension,
                adata_batch_list[i].obsp[adj_key])))
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
                adata_batch_list[i].obsp[adj_key],
                after_batch_connectivities_extension)))
        len_before_batch += adata_batch_list[i].shape[0]
    connectivities = sp.vstack(batch_connectivities)
    adata.obsp[adj_key] = connectivities
else:
    adata = ad.read_h5ad(
            f"{so_data_gold_folder_path}/{args.dataset}.h5ad")
    # Compute (separate) spatial neighborhood graphs
    sq.gr.spatial_neighbors(adata,
                            coord_type="generic",
                            spatial_key=args.spatial_key,
                            n_neighs=args.n_neighbors)
    # Make adjacency matrix symmetric
    adata.obsp[adj_key] = (
        adata.obsp[adj_key].maximum(
            adata.obsp[adj_key].T))
    
###############################################################################
### 3.2 Filter Genes ###
###############################################################################

print("\nFiltering for genes used in reference...")
adata = adata[:, genes]

###############################################################################
### 3.3 Add Gene Program Mask to Data ###
###############################################################################

# Add the gene program dictionary as binary masks to the adata for model 
# training. Use the same masks as for the reference model
adata.varm[gp_targets_mask_key] = reference_adata.varm[gp_targets_mask_key]
adata.varm[gp_sources_mask_key] = reference_adata.varm[gp_sources_mask_key]
adata.uns[args.gp_names_key] = reference_adata.uns[args.gp_names_key]
adata.uns[genes_idx_key] = reference_adata.uns[genes_idx_key]
adata.uns[target_genes_idx_key] = reference_adata.uns[target_genes_idx_key]
adata.uns[source_genes_idx_key] = reference_adata.uns[source_genes_idx_key]

###############################################################################
### 4. Model ###
###############################################################################

###############################################################################
## 4.1 Initialize, Train & Save Model ##
###############################################################################

print("\nTraining model...")
# Load model trained on reference data for transfer learning with query data
# Freeze all weights except for conditional weights
model = NicheCompass.load(
    dir_path=model_folder_path + f"/{args.reference_model_label}/" \
             f"{args.load_timestamp}",
    adata=adata,
    adata_file_name=f"{args.dataset}_{args.reference_model_label}.h5ad",
    gp_names_key=args.gp_names_key,
    unfreeze_all_weights=False,
    unfreeze_cond_embed_weights=True)

# Train model
model.train(n_epochs=args.n_epochs,
            n_epochs_all_gps=args.n_epochs_all_gps,
            n_epochs_no_cond_contrastive=args.n_epochs_no_cond_contrastive,
            lr=args.lr,
            lambda_edge_recon=args.lambda_edge_recon,
            lambda_gene_expr_recon=args.lambda_gene_expr_recon,
            lambda_cond_contrastive=args.lambda_cond_contrastive,
            contrastive_logits_pos_ratio=args.contrastive_logits_pos_ratio,
            contrastive_logits_neg_ratio=args.contrastive_logits_neg_ratio,
            lambda_group_lasso=args.lambda_group_lasso,
            lambda_l1_masked=args.lambda_l1_masked,
            edge_batch_size=args.edge_batch_size,
            node_batch_size=args.node_batch_size,
            mlflow_experiment_id=mlflow_experiment_id,
            verbose=True)

print("\nSaving query model...")
# Save trained model
model.save(
    dir_path=model_folder_path + f"/{args.query_model_label}/" \
             f"{args.load_timestamp}",
    overwrite=True,
    save_adata=True,
    adata_file_name=f"{args.dataset}_{args.query_model_label}.h5ad")

###############################################################################
## 4.2 Integrate AnnData, Add to Model & Save Model ##
###############################################################################

print("\nIntegrating reference and query adata...")
adata_batch_list = [reference_adata, model.adata]
reference_query_adata = ad.concat(adata_batch_list, join="inner")

# Combine spatial neighborhood graphs as disconnected components
batch_connectivities = []
len_before_batch = 0
for i in range(len(adata_batch_list)):
    if i == 0: # first batch
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[0].shape[0],
            (reference_query_adata.shape[0] -
            adata_batch_list[0].shape[0])))
        batch_connectivities.append(sp.hstack(
            (adata_batch_list[0].obsp[adj_key],
            after_batch_connectivities_extension)))
    elif i == (len(adata_batch_list) - 1): # last batch
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (reference_query_adata.shape[0] -
            adata_batch_list[i].shape[0])))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key])))
    else: # middle batches
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0], len_before_batch))
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (reference_query_adata.shape[0] -
            adata_batch_list[i].shape[0] -
            len_before_batch)))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key],
            after_batch_connectivities_extension)))
    len_before_batch += adata_batch_list[i].shape[0]
connectivities = sp.vstack(batch_connectivities)
reference_query_adata.obsp[adj_key] = connectivities

model.adata = reference_query_adata

model.adata.varm[gp_targets_mask_key] = reference_adata.varm[gp_targets_mask_key]
model.adata.varm[gp_sources_mask_key] = reference_adata.varm[gp_sources_mask_key]
model.adata.uns[args.gp_names_key] = reference_adata.uns[args.gp_names_key]
model.adata.uns[genes_idx_key] = reference_adata.uns[genes_idx_key]
model.adata.uns[target_genes_idx_key] = reference_adata.uns[target_genes_idx_key]
model.adata.uns[source_genes_idx_key] = reference_adata.uns[source_genes_idx_key]

print("\nComputing reference query latent embedding...")
model.adata.obsm[latent_key], _ = model.get_latent_representation(
   adata=model.adata,
   counts_key=counts_key,
   adj_key=adj_key,
   condition_key=condition_key,
   only_active_gps=True,
   return_mu_std=True,
   node_batch_size=model.node_batch_size_)

print("\nComputing neighbor graph...")
# Use latent representation for UMAP generation
sc.pp.neighbors(model.adata,
                use_rep=latent_key,
                key_added=latent_key)

print("\nComputing UMAP embedding...")
sc.tl.umap(model.adata,
           neighbors_key=latent_key)

# Store adata to disk
model.adata.write(
    f"{so_data_results_folder_path}/"
    f"{args.dataset}_nichecompass_{args.reference_query_model_label}.h5ad")

print("\nSaving reference query model...")
# Save trained model
model.save(
    dir_path=model_folder_path + f"/{args.reference_query_model_label}/" \
             f"{args.load_timestamp}",
    overwrite=True,
    save_adata=True,
    adata_file_name=f"{args.dataset}_{args.reference_query_model_label}.h5ad")