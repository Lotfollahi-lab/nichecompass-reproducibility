#!/bin/bash
#SBATCH -J deeplinc_single_sample_method_benchmarking_24
#SBATCH -o ../slurm_jobs/logs/out_deeplinc_single_sample_method_benchmarking_24.txt
#SBATCH -e ../slurm_jobs/logs/err_deeplinc_single_sample_method_benchmarking_24.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=128G
#SBATCH --nice=10000

dataset="slideseqv2_mouse_hippocampus_subsample_1pct"

source $HOME/.bashrc
conda activate deeplinc
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking/deeplinc
bash train_deeplinc_models.sh "$dataset"