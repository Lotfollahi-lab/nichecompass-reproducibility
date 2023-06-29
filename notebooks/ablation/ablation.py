import sys
sys.path.append("../../utils")

print("1")

import pandas as pd

from ablation_utils import *

ablation_task = "loss_weights"
datasets = ["seqfish_mouse_organogenesis", "xenium_human_breast_cancer"]
cell_type_keys = ["celltype_mapped_refined", "cell_states"]
condition_keys = [None, None]
experiment_ids = [2, 3]
latent_key = "nichecompass_latent"
spatial_key = "spatial"
latent_knng_key = "nichecompass_latent_knng"
spatial_knng_key = "spatial_knng"
gp_names_key = "nichecompass_gp_names"

artifact_folder_path = f"../../artifacts"

summary_df = pd.read_csv(f"mlflow_summary_{ablation_task}_ablation.csv")

print("2")

# Compute metrics and add to summary df
metrics_df = pd.DataFrame()
for i, dataset in enumerate(datasets):
    # Get timestamps of ablation runs for specific dataset
    timestamps = summary_df[(summary_df["dataset"] == dataset) & (summary_df["val_auroc_score"].notnull())]["timestamp"].tolist()
    
    # Compute metrics for ablation runs models
    current_iteration_metrics_df = compute_metrics(
        artifact_folder_path=artifact_folder_path,
        dataset=dataset,
        task=ablation_task + "_ablation",
        timestamps=timestamps,
        cell_type_key=cell_type_keys[i],
        condition_key=condition_keys[i],
        spatial_knng_key=spatial_knng_key,
        latent_knng_key=latent_knng_key,
        spatial_key=spatial_key,
        latent_key=latent_key)
    metrics_df = pd.concat([metrics_df, current_iteration_metrics_df], axis=0)
summary_df = pd.merge(summary_df, metrics_df, on=["dataset", "timestamp"], how="left")
summary_df.to_csv("temp.csv")