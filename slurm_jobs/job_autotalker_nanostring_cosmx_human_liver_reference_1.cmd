#!/bin/bash
#SBATCH -J autotalker_nanostring_cosmx_human_liver_reference_1
#SBATCH -o /home/aih/sebastian.birk/workspace/projects/autotalker-repro-new/slurm_jobs/logs/out_autotalker_nanostring_cosmx_human_liver_reference_1.txt
#SBATCH -e /home/aih/sebastian.birk/workspace/projects/autotalker-repro-new/slurm_jobs/logs/err_autotalker_nanostring_cosmx_human_liver_reference_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate autotalker_hpc
cd /
cd /home/aih/sebastian.birk/workspace/projects/autotalker-repro-new/scripts
python train_autotalker_reference_model.py --dataset=nanostring_cosmx_human_liver --batch=sample1 --nichenet_max_n_target_genes_per_gp=20000 --include_mebocost_gps=True --n_epochs=20 --n_epochs_all_gps=10 --lambda_group_lasso=0. --lambda_l1_masked=0. --edge_batch_size=256 --node_batch_size=32
