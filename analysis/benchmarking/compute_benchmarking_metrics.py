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

import anndata as ad
import pandas as pd
import scanpy as sc

root_folder_path = os.path.dirname(
    os.path.dirname(
        os.path.dirname(
            os.path.abspath(__file__))))

import sys
sys.path.append(f"{root_folder_path}/analysis/benchmarking/SDMBench/SDMBench")

from SDMBench import sdmbench

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
    "--niche_type_key",
    type=none_or_value,
    default=None,
    help="Niche type key corresponding to the dataset.")
parser.add_argument(
    "--batch_key",
    type=none_or_value,
    default="nichecompass_latent",
    help="Batch key corresponding to the dataset.")
parser.add_argument(
    "--batches",
    nargs='+',
    type=none_or_value,
    default=None,
    help="Batches.")
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
parser.add_argument(
    "--include_sdmbench",
    action=argparse.BooleanOptionalAction,
    default=False,
    help="Indicator whether to compute SDMBench metrics.")

args = parser.parse_args()

###############################################################################
### 1.3 Configure Paths ###
###############################################################################

artifact_folder_path = f"{root_folder_path}/artifacts"
so_data_folder_path = f"{root_folder_path}/datasets/st_data"
so_data_gold_folder_path = f"{so_data_folder_path}/gold"
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
        if metric in ["basw", "bgc", "blisi"]:
            continue
        else:
            benchmark_dict_acc[metric] = [] 
    else:
        benchmark_dict_acc[metric] = []
if args.include_sdmbench:
    for metric in ["sdmbench_ari",
                   "sdmbench_nmi",
                   "sdmbench_chaos",
                   "sdmbench_pas",
                   "sdmbench_asw",
                   "sdmbench_hom",
                   "sdmbench_com"]:
        benchmark_dict_acc[metric] = [] 

print("Loading adata}.")
adata = sc.read_h5ad(f"{file_folder_path}/{args.file_name}")

if ("pcr" in args.metrics) & (args.batches is not None):
    adata_batch_list = []
    for batch in args.batches:
        print(f"\nProcessing batch {batch}...")
        print("Loading data...")
        adata_batch = ad.read_h5ad(
            f"{so_data_gold_folder_path}/{args.dataset}_{batch}.h5ad")
        adata_batch_list.append(adata_batch)
    adata_raw = ad.concat(adata_batch_list, join="inner")
    sc.tl.pca(adata_raw, use_highly_variable=False)
    pcr_X_pre = adata_raw.obsm["X_pca"]
    del(adata_raw)
else:
    pcr_X_pre = None

for run_number in args.run_numbers:
    benchmark_dict_acc["dataset"].append(args.dataset)
    benchmark_dict_acc["run_number"].append(run_number)
    benchmark_dict_acc["run_time"].append(
        adata.uns[f"{args.latent_key.split('_')[0]}_model_training_duration_run{run_number}"])
    
    if args.include_sdmbench:
        # Compute latent neighbor graph
        sc.pp.neighbors(adata,
                        use_rep=f"{args.latent_key}_run{run_number}",
                        key_added=f"{args.latent_key}_run{run_number}")
        
        # Compute Leiden clustering of latent space until 'target_n_niches' niches are obtained (to match ground truth number)
        target_n_niches = adata.obs[args.niche_type_key].nunique()
        
        latent_leiden_resolution = 0.3
        latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"
        counter = 0
        while True:
            sc.tl.leiden(adata=adata,
                         resolution=latent_leiden_resolution,
                         key_added="pred_niche_types",
                         neighbors_key=f"{args.latent_key}_run{run_number}")

            niche_counts = adata.obs["pred_niche_types"].value_counts()
            valid_niches = niche_counts[niche_counts >= 100].index
            n_niches = adata.obs[adata.obs["pred_niche_types"].isin(valid_niches)]["pred_niche_types"].nunique()
            print(f"Current number of niches: {n_niches}")
            print(f"Cluster counter: {counter}")
            if n_niches == target_n_niches:
                break
            elif n_niches < target_n_niches and counter < 30:
                latent_leiden_resolution += 0.1
            elif n_niches < target_n_niches and counter > 30 and counter < 60:
                latent_leiden_resolution += 0.01
            elif n_niches > target_n_niches and counter < 30:
                latent_leiden_resolution -= 0.1
            elif n_niches > target_n_niches and counter > 30 and counter < 60:
                latent_leiden_resolution -= 0.01
            elif counter > 60:
                break
            counter += 1
        
        benchmark_dict_acc["sdmbench_ari"].append(sdmbench.compute_ARI(
            adata,
            args.niche_type_key,
            "pred_niche_types"))
        benchmark_dict_acc["sdmbench_nmi"].append(sdmbench.compute_NMI(
            adata,
            args.niche_type_key,
            "pred_niche_types"))
        benchmark_dict_acc["sdmbench_chaos"].append(sdmbench.compute_CHAOS(
            adata,
            "pred_niche_types"))
        benchmark_dict_acc["sdmbench_pas"].append(sdmbench.compute_PAS(
            adata,
            "pred_niche_types",
            spatial_key="spatial"))
        benchmark_dict_acc["sdmbench_asw"].append(sdmbench.compute_ASW(
            adata,
            "pred_niche_types",
            spatial_key=args.spatial_key))
        benchmark_dict_acc["sdmbench_hom"].append(sdmbench.compute_HOM(
            adata,
            args.niche_type_key,
            "pred_niche_types"))
        benchmark_dict_acc["sdmbench_com"].append(sdmbench.compute_COM(
            adata,
            args.niche_type_key,
            "pred_niche_types"))        
        
    benchmark_dict = compute_benchmarking_metrics(
            adata=adata,
            metrics=args.metrics,
            cell_type_key=args.cell_type_key,
            batch_key=args.batch_key,
            spatial_key=args.spatial_key,
            latent_key=f"{args.latent_key}_run{run_number}",
            pcr_X_pre=pcr_X_pre,
            n_jobs=1,
            seed=0,
            mlflow_experiment_id=None)

    for key, value in benchmark_dict.items():
        benchmark_dict_acc[key].append(value)

benchmark_df = pd.DataFrame(benchmark_dict_acc)

# Store results in temp file
benchmark_df.to_csv(f"{file_folder_path}/{args.file_name.replace('.h5ad', '_metrics.csv')}", index=False)