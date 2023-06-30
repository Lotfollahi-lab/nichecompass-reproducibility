#!/usr/bin/env python
# coding: utf-8

###############################################################################
# NicheCompass Metrics Computation #
###############################################################################

###############################################################################
## 1. Setup ##
###############################################################################

###############################################################################
### 1.1 Import Libraries ###
###############################################################################

import sys
sys.path.append("../../utils")
sys.path.append("../utils")

import argparse
import os

import pandas as pd

from ablation_utils import *

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

parser.add_argument(
    "--ablation_task",
    type=str,
    default="",
    help="Ablation task for which metrics will be computed.")
parser.add_argument(
    "--datasets",
    nargs="+",
    default="",
    help="Datasets for which metrics will be computed.")
parser.add_argument(
    "--cell_type_keys",
    nargs="+",
    default="",
    help="Cell type keys corresponding to the datasets.")
parser.add_argument(
    "--condition_keys",
    type=none_or_value,
    nargs="+",
    default="",
    help="Condition keys corresponding to the datasets.")
parser.add_argument(
    "--experiment_ids",
    nargs='+',
    type=none_or_int,
    default=None,
    help="Experiment IDs corresponding to the datasets")
parser.add_argument(
    "--latent_key",
    type=str,
    default="nichecompass_latent",
    help="s. NicheCompass metrics function signature.")
parser.add_argument(
    "--spatial_key",
    type=str,
    default="spatial",
    help="s. NicheCompass metrics function signature.")
parser.add_argument(
    "--latent_knng_key",
    type=str,
    default="nichecompass_latent_knng",
    help="s. NicheCompass metrics function signature.")
parser.add_argument(
    "--spatial_knng_key",
    type=str,
    default="spatial_knng",
    help="s. NicheCompass metrics function signature.")

args = parser.parse_args()

###############################################################################
### 1.3 Configure Paths ###
###############################################################################

root_folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
artifact_folder_path = f"{root_folder_path}/artifacts"
miscellaneous_folder_path = f"{artifact_folder_path}/miscellaneous"

###############################################################################
## 2. Compute Metrics ##
###############################################################################

summary_df = pd.read_csv(f"{miscellaneous_folder_path}/mlflow_summary_{args.ablation_task}_ablation.csv")

# Compute metrics and add to summary df
metrics_df = pd.DataFrame()
for i, dataset in enumerate(args.datasets):
    # Get timestamps of ablation runs for specific dataset
    timestamps = summary_df[(summary_df["dataset"] == dataset) & (summary_df["val_auroc_score"].notnull())]["timestamp"].tolist()
    
    # Compute metrics for ablation runs models
    current_iteration_metrics_df = compute_metrics(
        artifact_folder_path=artifact_folder_path,
        dataset=dataset,
        task=args.ablation_task + "_ablation",
        timestamps=timestamps,
        cell_type_key=args.cell_type_keys[i],
        condition_key=args.condition_keys[i],
        spatial_knng_key=args.spatial_knng_key,
        latent_knng_key=args.latent_knng_key,
        spatial_key=args.spatial_key,
        latent_key=args.latent_key)
    metrics_df = pd.concat([metrics_df, current_iteration_metrics_df], axis=0)
    metrics_df.to_csv(f"{miscellaneous_folder_path}/temp.csv")
summary_df = pd.merge(summary_df, metrics_df, on=["dataset", "timestamp"], how="left")
summary_df.to_csv(f"{miscellaneous_folder_path}/mlflow_summary_with_metrics_{ablation_task}_ablation.csv")