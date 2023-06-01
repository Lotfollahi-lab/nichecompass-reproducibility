#!/usr/bin/env python
# coding: utf-8

# # GraphST

# - **Creator**: Sebastian Birk (<sebastian.birk@helmholtz-munich.de>).
# - **Affiliation:** Helmholtz Munich, Institute of Computational Biology (ICB), Talavera-López Lab
# - **Date of Creation:** 11.01.2023
# - **Date of Last Modification:** 31.05.2023

# - The GraphST source code is available at https://github.com/JinmiaoChenLab/GraphST.
# - The corresponding preprint is "Long, Y. et al. DeepST: A versatile graph contrastive learning framework for spatially informed clustering, integration, and deconvolution of spatial transcriptomics. Preprint at https://doi.org/10.1101/2022.08.02.502407".
# - The workflow of this notebook follows the tutorial from https://deepst-tutorials.readthedocs.io/en/latest/Tutorial%204_Horizontal%20Integration.html.
# - The authors use raw counts as input to GraphST (stored in adata.X). Therefore, we also use raw counts.
# - To define the spatial neighborhood graph, the original GraphST paper uses the 3 nearest neighbors of a cell as neighbors and the union of all neighbors is used as final spatial neighborhood graph (the adjacency matrix is made symmetric). We use the same method but vary the number of neighbors between 4, 8, 12, 16 and 20.

# ## 1. Setup

# ### 1.1 Import Libraries

# In[1]:


import gc
import os
import time
from datetime import datetime

import anndata as ad
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import squidpy as sq
import torch
from GraphST import GraphST
from sklearn import metrics


# ### 1.2 Define Parameters

# In[2]:


model_name = "graphst"
latent_key = f"{model_name}_latent"
mapping_entity_key = "reference"
condition_key = "batch"
counts_key = "counts"
spatial_key = "spatial"
adj_key = "spatial_connectivities"
leiden_resolution = 0.3 # used for Leiden clustering of latent space
random_seed = 0 # used for Leiden clustering


# ### 1.3 Run Notebook Setup

# In[3]:


sc.set_figure_params(figsize=(6, 6))


# In[4]:


# Get time of notebook execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")


# In[5]:


# Run device. By default, the package is implemented on 'cpu'. It is recommended to use GPU.
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


# ### 1.4 Configure Paths and Directories

# In[6]:


srt_data_gold_folder_path = "../../datasets/srt_data/gold/"
srt_data_results_folder_path = "../../datasets/srt_data/results/" 
figure_folder_path = f"../../figures"

# Create required directories
os.makedirs(srt_data_gold_folder_path, exist_ok=True)
os.makedirs(srt_data_results_folder_path, exist_ok=True)


# ## 2. GraphST Model

# ### 2.1 Define Training Function

# In[7]:


def train_graphst_models(dataset,
                         reference_batches,
                         cell_type_key,
                         adata_new=None,
                         n_start_run=1,
                         n_end_run=10,
                         n_neighbor_list=[4, 4, 8, 8, 12, 12, 16, 16, 20, 20],
                         plot_latent_umaps: bool=False):    
    # Create new adata to store results from training runs in storage-efficient way
    if adata_new is None:  
        adata_batch_list = []
        if reference_batches is not None:
            for batch in reference_batches:
                adata_batch = ad.read_h5ad(
                    f"{srt_data_gold_folder_path}/{dataset}_{batch}.h5ad")
                adata_batch.obs[mapping_entity_key] = "reference"
                adata_batch_list.append(adata_batch)
            adata_original = ad.concat(adata_batch_list, join="inner")
        else:
            adata_original = ad.read_h5ad(f"{srt_data_gold_folder_path}/{dataset}.h5ad")

        adata_new = sc.AnnData(sp.csr_matrix(
            (adata_original.shape[0], adata_original.shape[1]),
            dtype=np.float32))
        adata_new.var_names = adata_original.var_names
        adata_new.obs_names = adata_original.obs_names
        adata_new.obs["cell_type"] = adata_original.obs[cell_type_key].values
        adata_new.obsm["spatial"] = adata_original.obsm["spatial"]
        adata_new.obs[condition_key] = adata_original.obs[condition_key]
        adata_new.obs[mapping_entity_key] = adata_original.obs[mapping_entity_key] 
        del(adata_original)

    model_seeds = list(range(10))
    for run_number, n_neighbors in zip(np.arange(n_start_run, n_end_run+1), n_neighbor_list):
        # Load data
        adata_batch_list = []
        if reference_batches is not None:
            for batch in reference_batches:
                print(f"Processing batch {batch}...")
                print("Loading data...")
                adata_batch = ad.read_h5ad(
                    f"{srt_data_gold_folder_path}/{dataset}_{batch}.h5ad")
                adata_batch.obs[mapping_entity_key] = "reference"
                print("Computing spatial neighborhood graph...\n")
                # Compute (separate) spatial neighborhood graphs
                sq.gr.spatial_neighbors(adata_batch,
                                        coord_type="generic",
                                        spatial_key=spatial_key,
                                        n_neighs=n_neighbors)
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
            adata = ad.read_h5ad(f"{srt_data_gold_folder_path}/{dataset}.h5ad")
            # Compute (separate) spatial neighborhood graphs
            sq.gr.spatial_neighbors(adata,
                                    coord_type="generic",
                                    spatial_key=spatial_key,
                                    n_neighs=n_neighbors)
            # Make adjacency matrix symmetric
            adata.obsp[adj_key] = (
                adata.obsp[adj_key].maximum(
                    adata.obsp[adj_key].T))
        
        # Store raw counts in adata.X
        adata.X = adata.layers["counts"]
        if "log1p" in adata.uns:
            del(adata.uns["log1p"])

        # Supply precomputed spatial neighborhood graph
        adata.obsm["graph_neigh"] = adata.obsp[adj_key]
        
        # Make adjacency matrix symmetric
        adata.obsm["adj"] = adata.obsm["graph_neigh"].maximum(
            adata.obsm["graph_neigh"].T)

        start_time = time.time()
        
        # Define model
        model = GraphST.GraphST(adata,
                                device=device,
                                random_seed=model_seeds[run_number-1])

        # Train model
        adata = model.train()
        
        # Measure time for model training
        end_time = time.time()
        elapsed_time = end_time - start_time
        hours, rem = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(rem, 60)
        print(f"Duration of model training in run {run_number}: "
              f"{int(hours)} hours, {int(minutes)} minutes and {int(seconds)} seconds.")
        adata_new.uns[f"{model_name}_model_training_duration_run{run_number}"] = (
            elapsed_time)

        if plot_latent_umaps:
            # Configure figure folder path
            dataset_figure_folder_path = f"{figure_folder_path}/{dataset}/sample_integration_method_benchmarking/" \
                                         f"{model_name}/{current_timestamp}"
            os.makedirs(dataset_figure_folder_path, exist_ok=True)
    
            # Use GraphST latent space for UMAP generation
            sc.pp.neighbors(adata,
                            use_rep="emb",
                            n_neighbors=n_neighbors)
            sc.tl.umap(adata)
            fig = sc.pl.umap(adata,
                             color=[cell_type_key],
                             title="Latent Space with Cell Types: GraphST",
                             return_fig=True)
            fig.savefig(f"{dataset_figure_folder_path}/latent_{model_name}"
                        f"_cell_types_run{run_number}.png",
                        bbox_inches="tight")

            # Compute latent Leiden clustering
            sc.tl.leiden(adata=adata,
                         resolution=leiden_resolution,
                         random_state=random_seed,
                         key_added=f"latent_graphst_leiden_{str(leiden_resolution)}")

            # Create subplot of latent Leiden cluster annotations in physical and latent space
            fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 12))
            title = fig.suptitle(t="Latent and Physical Space with Leiden Clusters: GraphST")
            sc.pl.umap(adata=adata,
                       color=[f"latent_{model_name}_leiden_{str(leiden_resolution)}"],
                       title=f"Latent Space with Leiden Clusters",
                       ax=axs[0],
                       show=False)
            sq.pl.spatial_scatter(adata=adata,
                                  color=[f"latent_{model_name}_leiden_{str(leiden_resolution)}"],
                                  title=f"Physical Space with Leiden Clusters",
                                  shape=None,
                                  ax=axs[1])

            # Create and position shared legend
            handles, labels = axs[0].get_legend_handles_labels()
            lgd = fig.legend(handles, labels, bbox_to_anchor=(1.25, 0.9185))
            axs[0].get_legend().remove()
            axs[1].get_legend().remove()

            # Adjust, save and display plot
            plt.subplots_adjust(wspace=0, hspace=0.2)
            fig.savefig(f"{dataset_figure_folder_path}/latent_physical_comparison_"
                        f"{model_name}_run{run_number}.png",
                        bbox_extra_artists=(lgd, title),
                        bbox_inches="tight")
            plt.show()

        # Store latent representation
        adata_new.obsm[latent_key + f"_run{run_number}"] = adata.obsm["emb"]
        
        # Use latent representation for UMAP generation
        sc.pp.neighbors(adata_new,
                        use_rep=f"{latent_key}_run{run_number}",
                        key_added=f"{latent_key}_run{run_number}")
        sc.tl.umap(adata_new,
                   neighbors_key=f"{latent_key}_run{run_number}")
        adata_new.obsm[f"{latent_key}_run{run_number}_X_umap"] = adata_new.obsm["X_umap"]
        del(adata_new.obsm["X_umap"])

        # Store intermediate adata to disk
        adata_new.write(f"{srt_data_results_folder_path}/sample_integration_method_benchmarking/"
                        f"{dataset}_{model_name}_sample_integration_method_benchmarking.h5ad")
        
        # Free memory
        del(adata)
        del(model)
        gc.collect()

    # Store final adata to disk
    adata_new.write(f"{srt_data_results_folder_path}/sample_integration_method_benchmarking/"
                    f"{dataset}_{model_name}_sample_integration_method_benchmarking.h5ad") 


# ### 2.2 Train Models on Benchmarking Datasets

# In[ ]:


train_graphst_models(dataset="seqfish_mouse_organogenesis",
                     reference_batches=[f"batch{i}" for i in range(1,7)],
                     cell_type_key="celltype_mapped_refined",
                     adata_new=None,
                     n_start_run=1,
                     n_end_run=10,
                     n_neighbor_list=[4, 4, 8, 8, 12, 12, 16, 16, 20, 20])


# In[ ]:


for subsample_pct in [1, 5, 10, 25, 50]:
    train_graphst_models(dataset=f"seqfish_mouse_organogenesis_subsample_{subsample_pct}pct",
                         reference_batches=[f"batch{i}" for i in range(1,7)],
                         cell_type_key="celltype_mapped_refined",
                         adata_new=None,
                         n_start_run=1,
                         n_end_run=10,
                         n_neighbor_list=[4, 4, 8, 8, 12, 12, 16, 16, 20, 20])


# In[ ]:


# This leads to memory overflow (NVIDIA GeForce RTX 3090 24GB GPU)
train_graphst_models(dataset="starmap_plus_mouse_cns",
                     reference_batches=[f"batch{i}" for i in range(1,21)],
                     cell_type_key="Main_molecular_cell_type",
                     adata_new=None,
                     n_start_run=1,
                     n_end_run=10,
                     n_neighbor_list=[4, 4, 8, 8, 12, 12, 16, 16, 20, 20])


# In[ ]:


# This leads to memory overflow starting at 5pct (NVIDIA GeForce RTX 3090 24GB GPU)
for subsample_pct in [1, 5, 10, 25, 50]:
    train_graphst_models(dataset=f"starmap_plus_mouse_cns_subsample_{subsample_pct}pct",
                         reference_batches=[f"batch{i}" for i in range(1,21)],
                         cell_type_key="Main_molecular_cell_type",
                         adata_new=None,
                         n_start_run=1,
                         n_end_run=10,
                         n_neighbor_list=[4, 4, 8, 8, 12, 12, 16, 16, 20, 20])


# In[ ]:




