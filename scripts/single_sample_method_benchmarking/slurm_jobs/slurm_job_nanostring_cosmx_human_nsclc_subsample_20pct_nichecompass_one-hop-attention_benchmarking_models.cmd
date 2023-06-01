#!/bin/bash
#SBATCH -J nanostring_cosmx_human_nsclc_subsample_20pct_nichecompass_one-hop-attention_benchmarking_models
#SBATCH -o /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking/slurm_jobs/logs/out_nanostring_cosmx_human_nsclc_subsample_20pct_nichecompass_one-hop-attention_benchmarking_models.txt
#SBATCH -e /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking/slurm_jobs/logs/err_nanostring_cosmx_human_nsclc_subsample_20pct_nichecompass_one-hop-attention_benchmarking_models.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking
bash train_nanostring_cosmx_human_nsclc_subsample_20pct_nichecompass_one-hop-attention_benchmarking_models.sh
