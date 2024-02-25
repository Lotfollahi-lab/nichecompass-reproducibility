#!/usr/bin/env python
# coding: utf-8

###############################################################################
# NicheCompass Reference Model Training #
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
from datetime import datetime

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from nichecompass.models import NicheCompass

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

parser.add_argument(
    "--dataset",
    type=str,
    help="Input dataset name. The adata file name has to be f'{dataset}.h5ad'.")
parser.add_argument(
    "--model_label",
    type=str,
    default="reference",
    help="Label of the model under which it was saved.")
parser.add_argument(
    "--load_timestamp",
    type=str,
    help="Timestamp of model to be loaded.")
parser.add_argument(
    "--include_atac_modality",
    nargs='+',
    type=none_or_bool,
    default=None,
    help="s. NicheCompass class signature")
parser.add_argument(
    "--gp_names_key",
    type=str,
    default="nichecompass_gp_names",
    help="s. NicheCompass class signature")
parser.add_argument(
    "--node_batch_size",
    type=int,
    default=32768,
    help="s. NicheCompass train method signature")
parser.add_argument(
    "--compute_latent",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether to compute the latent representations.")
parser.add_argument(
    "--compute_pca",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether to compute PCA.")
parser.add_argument(
    "--pca_key",
    type=str,
    default="nichecompass_pca",
    help="Key in `adata.obsm` under which the PCA representations are stored.")
parser.add_argument(
    "--use_pca_knn_graph",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="Indicator whether to compute knn graph from PCA representations.")
parser.add_argument(
    "--compute_knn_graph",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether to compute a knn graph.")
parser.add_argument(
    "--compute_umap",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether to compute UMAP representations.")
parser.add_argument(
    "--compute_leiden",
    action=argparse.BooleanOptionalAction,
    default=True,
    help="Indicator whether to compute Leiden clustering.")
parser.add_argument(
    "--latent_leiden_resolution",
    type=float,
    default=0.2,
    help="Resolution used for Leiden clustering.")

args = parser.parse_args()

if args.include_atac_modality:
    save_adata_atac = True
else:
    save_adata_atac = False
    
latent_cluster_key = f"latent_leiden_{str(args.latent_leiden_resolution)}"

# Get time of script execution
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")

print(f"Run timestamp: {args.load_timestamp}.")
print("Script arguments:")
print(sys.argv)

###############################################################################
### 1.3 Configure Paths and Create Directories ###
###############################################################################

root_folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
artifacts_folder_path = f"{root_folder_path}/artifacts"
model_folder_path = f"{artifacts_folder_path}/{args.dataset}/models/" \
                    f"{args.model_label}/{args.load_timestamp}"

###############################################################################
### 2. Load Model & Perform Postprocessing ###
###############################################################################

print("\nLoading model...")
if args.include_atac_modality:
    model = NicheCompass.load(dir_path=model_folder_path,
                              adata=None,
                              adata_file_name=f"{args.dataset}_{args.model_label}.h5ad",
                              adata_atac=None,
                              adata_atac_file_name=f"{args.dataset}_{args.model_label}_atac.h5ad",
                              gp_names_key=args.gp_names_key)
else:
    model = NicheCompass.load(dir_path=model_folder_path,
                              adata=None,
                              adata_file_name=f"{args.dataset}_{args.model_label}.h5ad",
                              gp_names_key=args.gp_names_key)
    
if args.compute_latent:
    print("\nComputing latent representations...")
    model.adata.obsm[model.latent_key_] = model.get_latent_representation(
       adata=model.adata,
       counts_key=model.counts_key_,
       adj_key=model.adj_key_,
       cat_covariates_keys=model.cat_covariates_keys_,
       only_active_gps=True,
       return_mu_std=False,
       node_batch_size=args.node_batch_size,
       dtype=np.float16)
    print("\nSkipping latent representations computation...")
    
if args.compute_pca:
    print("\nComputing PCA representations...")
    model.adata.obsm[args.pca_key] = sc.tl.pca(model.adata.obsm[model.latent_key_])
else:
    print("\Skipping PCA representations computation...")
    
if args.compute_knn_graph:
    if args.use_pca_knn_graph:
        neighbor_rep_key = args.pca_key
    else:
        neighbor_rep_key = model.latent_key_
    print("\nComputing neighbor graph on PCA representations...")
    sc.pp.neighbors(model.adata,
                    use_rep=neighbor_rep_key,
                    key_added=model.latent_key_)
    print("\nSkipping neighbor graph computation...")

if args.compute_umap:
    print("\nComputing UMAP embedding...")
    sc.tl.umap(model.adata,
               neighbors_key=model.latent_key_)
    print("\nSkipping UMAP embedding computation...")
    
if args.compute_leiden:
    print("\nComputing Leiden clustering...")
    sc.tl.leiden(adata=model.adata,
                 resolution=args.latent_leiden_resolution,
                 key_added=latent_cluster_key,
                 neighbors_key=model.latent_key_)
    print("\nSkipping Leiden clustering computation...")

print("\nSaving model...")
model.save(dir_path=model_folder_path,
           overwrite=True,
           save_adata=True,
           adata_file_name=f"{args.dataset}_{args.model_label}.h5ad",
           save_adata_atac=save_adata_atac,
           adata_atac_file_name=f"{args.dataset}_{args.model_label}_atac.h5ad")
