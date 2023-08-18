#!/bin/bash
#SBATCH -J deeplinc_single_sample_method_benchmarking_7
#SBATCH -o ../slurm_jobs/logs/out_deeplinc_single_sample_method_benchmarking_7.txt
#SBATCH -e ../slurm_jobs/logs/err_deeplinc_single_sample_method_benchmarking_7.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=128G
#SBATCH --nice=10000

dataset="nanostring_cosmx_human_nsclc_batch5"

source $HOME/.bashrc
conda activate deeplinc
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking/deeplinc
bash train_deeplinc_models.sh "$dataset"