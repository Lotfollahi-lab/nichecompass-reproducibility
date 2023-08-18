#!/usr/bin/env python
# coding: utf-8

###############################################################################
# NicheCompass Benchmarking Metrics Computation #
###############################################################################

###############################################################################
## 1. Setup ##
###############################################################################

###############################################################################
### 1.1 Import Libraries ###
###############################################################################

import argparse
import os

# has to be before scanpy
from nichecompass.benchmarking import compute_benchmarking_metrics

import pandas as pd
import scanpy as sc

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
    "--dataset",
    type=str,
    default="",
    help="Dataset for which metrics will be computed.")
parser.add_argument(
    "--task",
    type=str,
    default="",
    help="Task for which metrics will be computed.")
parser.add_argument(
    "--file_name",
    type=str,
    default="",
    help="File name of adata.")
parser.add_argument(
    "--cell_type_key",
    type=str,
    default="cell_type",
    help="Cell type key corresponding to the dataset.")
parser.add_argument(
    "--batch_key",
    type=none_or_value,
    default="nichecompass_latent",
    help="Batch key corresponding to the dataset.")
parser.add_argument(
    "--spatial_key",
    type=str,
    default="spatial",
    help="s. NicheCompass metrics function signature.")
parser.add_argument(
    "--latent_key",
    type=str,
    default="nichecompass_latent",
    help="s. NicheCompass metrics function signature.")
parser.add_argument(
    "--run_numbers",
    nargs="+",
    default=[1, 2, 3, 4, 5, 6, 7, 8],
    help="Run numbers for which metrics should be computed.")
parser.add_argument(
    "--metrics",
    nargs="+",
    default=["nasw"],
    help="Metrics to be computed.")

args = parser.parse_args()

###############################################################################
### 1.3 Configure Paths ###
###############################################################################

root_folder_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
artifact_folder_path = f"{root_folder_path}/artifacts"
file_folder_path = f"{artifact_folder_path}/{args.task}"

###############################################################################
## 2. Compute Metrics ##
###############################################################################

# Initialize metrics dicts
benchmark_dict_acc = {"dataset": [],
                      "run_number": [],
                      "run_time": []}
for metric in args.metrics:
    if args.batch_key is None:
        if metric in ["basw", "bgc", "bilisi"]:
            continue
        else:
            benchmark_dict_acc[metric] = [] 
    else:
        benchmark_dict_acc[metric] = []

print("Loading adata}.")
adata = sc.read_h5ad(f"{file_folder_path}/{args.file_name}")

for run_number in args.run_numbers:
    benchmark_dict_acc["dataset"].append(args.dataset)
    benchmark_dict_acc["run_number"].append(run_number)
    benchmark_dict_acc["run_time"].append(
        adata.uns[f"{args.latent_key.split('_')[0]}_model_training_duration_"
                  f"run{run_number}"])
    
    benchmark_dict = compute_benchmarking_metrics(
            adata=adata,
            metrics=args.metrics,
            cell_type_key=args.cell_type_key,
            batch_key=args.batch_key,
            spatial_key=args.spatial_key,
            latent_key=f"{args.latent_key}_run{run_number}",
            n_jobs=1,
            seed=0,
            mlflow_experiment_id=None)

    for key, value in benchmark_dict.items():
        benchmark_dict_acc[key].append(value)

benchmark_df = pd.DataFrame(benchmark_dict_acc)

# Store results in temp file
benchmark_df.to_csv(f"{file_folder_path}/{args.file_name.replace('.h5ad', '_metrics.csv')}", index=False)