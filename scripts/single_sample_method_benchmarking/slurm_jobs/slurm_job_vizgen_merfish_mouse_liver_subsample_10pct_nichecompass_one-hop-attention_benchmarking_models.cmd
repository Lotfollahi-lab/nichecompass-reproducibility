#!/bin/bash
#SBATCH -J vizgen_merfish_mouse_liver_subsample_10pct_nichecompass_one-hop-attention_benchmarking_models
#SBATCH -o /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking/slurm_jobs/logs/out_vizgen_merfish_mouse_liver_subsample_10pct_nichecompass_one-hop-attention_benchmarking_models.txt
#SBATCH -e /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking/slurm_jobs/logs/err_vizgen_merfish_mouse_liver_subsample_10pct_nichecompass_one-hop-attention_benchmarking_models.txt
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
bash train_vizgen_merfish_mouse_liver_subsample_10pct_nichecompass_one-hop-attention_benchmarking_models.sh
