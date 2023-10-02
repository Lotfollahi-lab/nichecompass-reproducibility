# Imports =================
print('imports')
import sys
sys.path.append("../utils")

import argparse
import gc
import os
import random
import shutil
import warnings
from datetime import datetime
from matplotlib import rcParams

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import mlflow
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats
import seaborn as sns
import squidpy as sq
import torch
from matplotlib import gridspec
from matplotlib.pyplot import rc_context
from sklearn.preprocessing import MinMaxScaler
from pywaffle import Waffle

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                aggregate_obsp_matrix_per_cell_type,
                                create_cell_type_chord_plot_from_df,
                                create_new_color_dict,
                                generate_enriched_gp_info_plots)

from analysis_utils import (add_cell_type_latent_cluster_emphasis,
                            add_sub_cell_type,
                            compute_cell_type_latent_clusters,
                            generate_gp_info_plots,
                            plot_physical_latent_for_cell_types,
                            plot_cell_type_latent_clusters,
                            plot_latent,
                            plot_category_in_latent_and_physical_space,
                            sankey,
                            store_top_gps_summary)

# Parameters =================
print('define parameters')

dataset = "nanostring_cosmx_human_nsclc"

## AnnData keys
adj_key = "spatial_connectivities"
spatial_key = "spatial"
sub_cell_type_key = "cell_type_original"
nicke_key = "niche"
gp_names_key = "nichecompass_gp_names"
active_gp_names_key = "nichecompass_active_gp_names"
latent_key = "nichecompass_latent"
mapping_entity_key = "mapping_entity"
## Analysis
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"
## Others
random_seed = 0
multimodal = False
log_norm_omics_features = False
cell_type_groups = []
latent_groups = []
load_timestamp = "03092023_001459_4"
model_label = "reference_query_mapping"
latent_leiden_resolution = 0.5
latent_cluster_spot_size = 0.03
dataset_str = "nanoString CosMx Human NSCLC"
condition_key = "batch"
sample_key = "batch"
spot_size = 30
cell_type_key = "cell_type"
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"

# Figures
sc.set_figure_params(figsize=(6, 6))
sns.set_style("whitegrid", {'axes.grid' : False})

# Define paths
figure_folder_path = f"../artifacts/{dataset}/figures/{model_label}/{load_timestamp}"
model_folder_path = f"../artifacts/{dataset}/models/{model_label}/{load_timestamp}"
result_folder_path = f"../artifacts/{dataset}/results/{model_label}/{load_timestamp}"
gp_data_folder_path = "../datasets/gp_data" # gene program data
srt_data_folder_path = "../datasets/srt_data" # spatially resolved transcriptomics data
srt_data_gold_folder_path = f"{srt_data_folder_path}/gold"

# Create required directories
os.makedirs(figure_folder_path, exist_ok=True)
os.makedirs(result_folder_path, exist_ok=True)

# Load =================
print('load model')

model = NicheCompass.load(dir_path=model_folder_path,
                          adata=None,
                          adata_file_name=f"{dataset}_{model_label}_postprocessed.h5ad",
                          gp_names_key=gp_names_key)

selected_cats = None
log_bayes_factor_thresh = 2.3 # 2.3 strong threshold; 4.6 decisive threshold (https://en.wikipedia.org/wiki/Bayes_factor)
title = f"NicheCompass Latent Cluster Enriched Gene Programs Log Bayes Factor {log_bayes_factor_thresh}"
save_fig = True
file_path = f"{figure_folder_path}/res_{latent_leiden_resolution}_" \
            f"latent_clusters_all_vs_rest_log_bayes_factor_" \
            f"{log_bayes_factor_thresh}_enriched_gps_heatmap.pdf"

# Run differential gp testing
print('differential testing')
# Imports =================
print('imports', flush=True)
import sys
sys.path.append("../utils")

import argparse
import gc
import os
import random
import shutil
import warnings
from datetime import datetime
from matplotlib import rcParams

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import mlflow
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats
import seaborn as sns
import squidpy as sq
import torch
from matplotlib import gridspec
from matplotlib.pyplot import rc_context
from sklearn.preprocessing import MinMaxScaler
from pywaffle import Waffle

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                aggregate_obsp_matrix_per_cell_type,
                                create_cell_type_chord_plot_from_df,
                                create_new_color_dict,
                                generate_enriched_gp_info_plots)

from analysis_utils import (add_cell_type_latent_cluster_emphasis,
                            add_sub_cell_type,
                            compute_cell_type_latent_clusters,
                            generate_gp_info_plots,
                            plot_physical_latent_for_cell_types,
                            plot_cell_type_latent_clusters,
                            plot_latent,
                            plot_category_in_latent_and_physical_space,
                            sankey,
                            store_top_gps_summary)

# Parameters =================
print('define parameters', flush=True)

dataset = "nanostring_cosmx_human_nsclc"

## AnnData keys
adj_key = "spatial_connectivities"
spatial_key = "spatial"
sub_cell_type_key = "cell_type_original"
nicke_key = "niche"
gp_names_key = "nichecompass_gp_names"
active_gp_names_key = "nichecompass_active_gp_names"
latent_key = "nichecompass_latent"
mapping_entity_key = "mapping_entity"
## Analysis
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"
## Others
random_seed = 0
multimodal = False
log_norm_omics_features = False
cell_type_groups = []
latent_groups = []
load_timestamp = "03092023_001459_4"
model_label = "reference_query_mapping"
latent_leiden_resolution = 0.5
latent_cluster_spot_size = 0.03
dataset_str = "nanoString CosMx Human NSCLC"
condition_key = "batch"
sample_key = "batch"
spot_size = 30
cell_type_key = "cell_type"
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"

# Figures
sc.set_figure_params(figsize=(6, 6))
sns.set_style("whitegrid", {'axes.grid' : False})

# Define paths
figure_folder_path = f"../artifacts/{dataset}/figures/{model_label}/{load_timestamp}"
model_folder_path = f"../artifacts/{dataset}/models/{model_label}/{load_timestamp}"
result_folder_path = f"../artifacts/{dataset}/results/{model_label}/{load_timestamp}"
gp_data_folder_path = "../datasets/gp_data" # gene program data
srt_data_folder_path = "../datasets/srt_data" # spatially resolved transcriptomics data
srt_data_gold_folder_path = f"{srt_data_folder_path}/gold"

# Create required directories
os.makedirs(figure_folder_path, exist_ok=True)
os.makedirs(result_folder_path, exist_ok=True)

# Load =================
print('load model', flush=True)

model = NicheCompass.load(dir_path=model_folder_path,
                          adata=None,
                          adata_file_name=f"{dataset}_{model_label}_postprocessed.h5ad",
                          gp_names_key=gp_names_key)

selected_cats = None
log_bayes_factor_thresh = 2.3 # 2.3 strong threshold; 4.6 decisive threshold (https://en.wikipedia.org/wiki/Bayes_factor)
title = f"NicheCompass Latent Cluster Enriched Gene Programs Log Bayes Factor {log_bayes_factor_thresh}"
save_fig = True
file_path = f"{figure_folder_path}/res_{latent_leiden_resolution}_" \
            f"latent_clusters_all_vs_rest_log_bayes_factor_" \
            f"{log_bayes_factor_thresh}_enriched_gps_heatmap.pdf"

# Run differential gp testing
print('differential testing')
gp_summary_df = model.get_gp_summary()

enriched_gps = model.run_differential_gp_tests(
    cat_key=latent_cluster_key,
    selected_cats=selected_cats,
    comparison_cats="rest",
    log_bayes_factor_thresh=log_bayes_factor_thresh)

print('save', flush=True)

model.adata.write(f"{model_folder_path}/{dataset}_{model_label}_postprocessed.h5ad")

save_file = True
file_path = f"{figure_folder_path}/res_{latent_leiden_resolution}_" \
            f"latent_clusters_all_vs_rest_log_bayes_factor_" \
            f"{log_bayes_factor_thresh}_enriched_gps_summary.csv"

gp_summary_cols = ["gp_name",
                   "n_source_genes",
                   "n_non_zero_source_genes",
                   "n_target_genes",
                   "n_non_zero_target_genes",
                   "gp_source_genes",
                   "gp_target_genes",
                   "gp_source_genes_weights",
                   "gp_target_genes_weights",
                   "gp_source_genes_importances",
                   "gp_target_genes_importances"]

# Get summary of decisively enriched gene programs
enriched_gp_summary_df = gp_summary_df[gp_summary_df["gp_name"].isin(enriched_gps)]
cat_dtype = pd.CategoricalDtype(categories=enriched_gps, ordered=True)
enriched_gp_summary_df["gp_name"] = enriched_gp_summary_df["gp_name"].astype(cat_dtype)
enriched_gp_summary_df = enriched_gp_summary_df.sort_values(by="gp_name")
enriched_gp_summary_df = enriched_gp_summary_df[gp_summary_cols]

if save_file:
    enriched_gp_summary_df.to_csv(f"{file_path}")
else:
    display(enriched_gp_summary_df)


# Plot heatmap of enriched gps
print('plot', flush=True)

ax = sc.pl.heatmap(model.adata,
                    enriched_gps,
                    show_gene_labels=True,
                    groupby=latent_cluster_key,
                    dendrogram=True,
                    swap_axes=True,
                    figsize=(model.adata.obs[latent_cluster_key].nunique() * 1.5,
                             len(enriched_gps) / 2),
                    save=save_fig)
if save_fig:
    shutil.move("figures/heatmap.pdf", file_path)
    os.rmdir("figures")
enriched_gps = model.run_differential_gp_tests(
    cat_key=latent_cluster_key,
    selected_cats=selected_cats,
    comparison_cats="rest",
    log_bayes_factor_thresh=log_bayes_factor_thresh)

print('save')

model.adata.write(f"{model_folder_path}/{dataset}_{model_label}_postprocessed.h5ad")

save_file = True
file_path = f"{figure_folder_path}/res_{latent_leiden_resolution}_" \
            f"latent_clusters_all_vs_rest_log_bayes_factor_" \
            f"{log_bayes_factor_thresh}_enriched_gps_summary.csv"

gp_summary_cols = ["gp_name",
                   "n_source_genes",
                   "n_non_zero_source_genes",
                   "n_target_genes",
                   "n_non_zero_target_genes",
                   "gp_source_genes",
                   "gp_target_genes",
                   "gp_source_genes_weights",
                   "gp_target_genes_weights",
                   "gp_source_genes_importances",
                   "gp_target_genes_importances"]

# Get summary of decisively enriched gene programs
enriched_gp_summary_df = gp_summary_df[gp_summary_df["gp_name"].isin(enriched_gps)]
cat_dtype = pd.CategoricalDtype(categories=enriched_gps, ordered=True)
enriched_gp_summary_df["gp_name"] = enriched_gp_summary_df["gp_name"].astype(cat_dtype)
enriched_gp_summary_df = enriched_gp_summary_df.sort_values(by="gp_name")
enriched_gp_summary_df = enriched_gp_summary_df[gp_summary_cols]

if save_file:
    enriched_gp_summary_df.to_csv(f"{file_path}")
else:
    display(enriched_gp_summary_df)


# Plot heatmap of enriched gps
print('plot')

ax = sc.pl.heatmap(model.adata,
                    enriched_gps,
                    show_gene_labels=True,
                    groupby=latent_cluster_key,
                    dendrogram=True,
                    swap_axes=True,
                    figsize=(model.adata.obs[latent_cluster_key].nunique() * 1.5,
                             len(enriched_gps) / 2),
                    save=save_fig)
if save_fig:
    shutil.move("figures/heatmap.pdf", file_path)
    os.rmdir("figures")
