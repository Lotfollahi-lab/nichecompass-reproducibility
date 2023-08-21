#!/usr/bin/env python
# coding: utf-8

# # expiMap

# - **Creator**: Sebastian Birk (<sebastian.birk@helmholtz-munich.de>).
# - **Affiliation:** Helmholtz Munich, Institute of Computational Biology (ICB), Talavera-López Lab
# - **Date of Creation:** 05.01.2023
# - **Date of Last Modification:** 19.08.2023

# - The expiMap source code is available at https://github.com/theislab/scarches.
# - The corresponding preprint is "Lotfollahi, M. et al. Biologically informed deep learning to infer gene program activity in single cells. bioRxiv 2022.02.05.479217 (2022) doi:10.1101/2022.02.05.479217".
# - The workflow of this notebook follows the tutorial from https://scarches.readthedocs.io/en/latest/expimap_surgery_pipeline_basic.html.
# - We use a modified version of the NicheCompass gene program mask with only target genes as the gene program mask for expimap. The reasons are that it is relevant for cell communication, to improve comparability and since the expiMap method did not work well on this dataset with the reactome gene program used in the above cited tutorial.
# - The authors use raw counts as input to expiMap. Therefore, we also use raw counts (stored in adata.X).

# ## 1. Setup

# ### 1.1 Import Libraries

# In[ ]:


import os
import time
from datetime import datetime

import gdown
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scarches as sca
import scipy.sparse as sp
import squidpy as sq
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                extract_gp_dict_from_mebocost_es_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                filter_and_combine_gp_dict_gps)


# ### 1.2 Define Parameters

# In[ ]:


model_name = "expimap"
latent_key = f"{model_name}_latent"
leiden_resolution = 0.5 # used for Leiden clustering of latent space
random_seed = 0 # used for Leiden clustering


# ### 1.3 Run Notebook Setup

# In[ ]:


sc.set_figure_params(figsize=(6, 6))


# In[ ]:


# Get time of notebook execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%d%m%Y_%H%M%S")


# ### 1.4 Configure Paths and Directories

# In[ ]:


data_folder_path = "../../datasets/srt_data/gold/"
benchmarking_folder_path = "../../artifacts/single_sample_method_benchmarking"
figure_folder_path = f"../../figures"
gp_data_folder_path = "../../datasets/gp_data" # gene program data
ga_data_folder_path = "../../datasets/ga_data" # gene annotation data

# Create required directories
os.makedirs(gp_data_folder_path, exist_ok=True)


# ## 2. expiMap Model

# ### 2.1 Prepare Gene Program Mask

# #### 2.1.1 Mouse

# In[ ]:


species = "mouse"

nichenet_lr_network_file_path = gp_data_folder_path + f"/nichenet_lr_network_v2_{species}.csv"
nichenet_ligand_target_matrix_file_path = gp_data_folder_path + f"/nichenet_ligand_target_matrix_v2_{species}.csv"
omnipath_lr_network_file_path = gp_data_folder_path + "/omnipath_lr_network.csv"
gene_orthologs_mapping_file_path = ga_data_folder_path + "/human_mouse_gene_orthologs.csv"

print("\nPreparing the gene program mask...")
# OmniPath gene programs
mouse_omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
    species=species,
    min_curation_effort=0,
    load_from_disk=True,
    save_to_disk=False,
    lr_network_file_path=omnipath_lr_network_file_path,
    gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
    plot_gp_gene_count_distributions=False)

# NicheNet gene programs
mouse_nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
    species=species,
    version="v2",
    keep_target_genes_ratio=1.0,
    max_n_target_genes_per_gp=250,
    load_from_disk=True,
    save_to_disk=False,
    lr_network_file_path=nichenet_lr_network_file_path,
    ligand_target_matrix_file_path=nichenet_ligand_target_matrix_file_path,
    gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
    plot_gp_gene_count_distributions=False)

# Combine gene programs into one dictionary
mouse_combined_gp_dict = dict(mouse_omnipath_gp_dict)
mouse_combined_gp_dict.update(mouse_nichenet_gp_dict)

mouse_mebocost_gp_dict = extract_gp_dict_from_mebocost_es_interactions(
    dir_path=f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps",
    species=species,
    plot_gp_gene_count_distributions=False)

mouse_combined_gp_dict.update(mouse_mebocost_gp_dict)
    
# Filter and combine gene programs
mouse_combined_new_gp_dict = filter_and_combine_gp_dict_gps(
    gp_dict=mouse_combined_gp_dict,
    gp_filter_mode="subset",
    combine_overlap_gps=True,
    overlap_thresh_source_genes=0.9,
    overlap_thresh_target_genes=0.9,
    overlap_thresh_genes=0.9,
    verbose=False)

print("Number of gene programs before filtering and combining: "
      f"{len(mouse_combined_new_gp_dict)}.")
print(f"Number of gene programs after filtering and combining: "
      f"{len(mouse_combined_new_gp_dict)}.")


# ### 2.2 Define Training Function

# In[ ]:


def train_expimap_models(dataset,
                         gp_dict,
                         cell_type_key,
                         adata_new=None,
                         n_start_run=1,
                         n_end_run=8,
                         n_neighbor_list=[4, 4, 8, 8, 12, 12, 16, 16],
                         plot_latent_umaps: bool=False):
    
    # Configure figure folder path
    dataset_figure_folder_path = f"{figure_folder_path}/{dataset}/method_benchmarking/expimap/{current_timestamp}"
    os.makedirs(dataset_figure_folder_path, exist_ok=True)
    
    # Create new adata to store results from training runs in storage-efficient way
    if adata_new is None:
        adata_original = sc.read_h5ad(data_folder_path + f"{dataset}.h5ad")
        adata_new = sc.AnnData(sp.csr_matrix(
            (adata_original.shape[0], adata_original.shape[1]),
            dtype=np.float32))
        adata_new.var_names = adata_original.var_names
        adata_new.obs_names = adata_original.obs_names
        adata_new.obs["cell_type"] = adata_original.obs[cell_type_key].values
        adata_new.obsm["spatial"] = adata_original.obsm["spatial"]
        del(adata_original)
    
    model_seeds = list(range(10))
    for run_number, n_neighbors in zip(np.arange(n_start_run, n_end_run+1), n_neighbor_list):
        # n_neighbors is here only used for the latent neighbor graph construction used for
        # UMAP generation and clustering as expiMap is not a spatial method
        
        # Load data
        adata = sc.read_h5ad(data_folder_path + f"{dataset}.h5ad")
        
        # Store raw counts in optimized format in adata.X
        adata.layers["counts"] = adata.layers["counts"].tocsr()
        adata.X = adata.layers["counts"]
        
        adata.obs["batch"] == "batch1"  
        
        # Add the gene program dictionary as binary masks to the adata for model training
        # Use only target genes from the NicheCompass gene program mask
        add_gps_from_gp_dict_to_adata(
            gp_dict=gp_dict,
            adata=adata,
            genes_uppercase=True,
            gp_targets_mask_key="I",
            gp_sources_mask_key="_",
            gp_names_key="terms",
            min_genes_per_gp=1,
            min_source_genes_per_gp=0,
            min_target_genes_per_gp=0,
            max_genes_per_gp=None,
            max_source_genes_per_gp=None,
            max_target_genes_per_gp=None)

        # Determine dimensionality of hidden encoder
        n_hidden_encoder = len(adata.uns["terms"])
        
        start_time = time.time()
        
        # Initialize model
        intr_cvae = sca.models.EXPIMAP(adata=adata,
                                       condition_key="batch",
                                       hidden_layer_sizes=[256, 256, 256],
                                       recon_loss="nb")

        # Train model
        early_stopping_kwargs = {
            "early_stopping_metric": "val_unweighted_loss",
            "threshold": 0,
            "patience": 50,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1}
        intr_cvae.train(
            n_epochs=400,
            alpha_epoch_anneal=100,
            alpha=0.7,
            alpha_kl=0.5,
            weight_decay=0.,
            early_stopping_kwargs=early_stopping_kwargs,
            use_early_stopping=True,
            monitor_only_val=False,
            seed=model_seeds[run_number-1])

        # Store latent representation
        adata.obsm[latent_key] = intr_cvae.get_latent(mean=False, only_active=True)
        
        # Measure time for model training
        end_time = time.time()
        elapsed_time = end_time - start_time
        hours, rem = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(rem, 60)
        print(f"Duration of model training in run {run_number}: {int(hours)} hours, {int(minutes)} minutes and {int(seconds)} seconds.")
        adata_new.uns[f"{model_name}_model_training_duration_run{run_number}"] = (
            elapsed_time)

        if plot_latent_umaps:
            # Use expiMap latent space for UMAP generation
            sc.pp.neighbors(adata,
                            use_rep=latent_key,
                            n_neighbors=n_neighbors)
            sc.tl.umap(adata)
            fig = sc.pl.umap(adata,
                             color=[cell_type_key],
                             title="Latent Space with Cell Types: expiMap",
                             return_fig=True)
            fig.savefig(f"{dataset_figure_folder_path}/latent_{model_name}"
                        f"_cell_types_run{run_number}.png",
                        bbox_inches="tight")

            # Compute latent Leiden clustering
            sc.tl.leiden(adata=adata,
                         resolution=leiden_resolution,
                         random_state=random_seed,
                         key_added=f"latent_{model_name}_leiden_{str(leiden_resolution)}")

            # Create subplot of latent Leiden cluster annotations in physical and latent space
            fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6, 12))
            title = fig.suptitle(t="Latent and Physical Space with Leiden Clusters: expiMap")
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
        adata_new.obsm[latent_key + f"_run{run_number}"] = adata.obsm[latent_key]

        # Store intermediate adata to disk
        adata_new.write(f"{benchmarking_folder_path}/{dataset}_{model_name}8.h5ad")

    # Store final adata to disk
    adata_new.write(f"{benchmarking_folder_path}/{dataset}_{model_name}8.h5ad")    


train_expimap_models(dataset="vizgen_merfish_mouse_liver_subsample_50pct",
                     gp_dict=mouse_combined_new_gp_dict,
                     cell_type_key="Cell_Type",
                     adata_new=None,
                     n_start_run=5,
                     n_end_run=8,
                     n_neighbor_list=[12, 12, 16, 16])
