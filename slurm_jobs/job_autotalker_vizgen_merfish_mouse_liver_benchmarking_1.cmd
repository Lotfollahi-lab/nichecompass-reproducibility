#!/bin/bash
#SBATCH -J autotalker_vizgen_merfish_mouse_liver_benchmarking_1
#SBATCH -o /home/aih/sebastian.birk/workspace/projects/autotalker-repro-new/slurm_jobs/logs/out_autotalker_vizgen_merfish_mouse_liver_benchmarking_1.txt
#SBATCH -e /home/aih/sebastian.birk/workspace/projects/autotalker-repro-new/slurm_jobs/logs/err_autotalker_vizgen_merfish_mouse_liver_benchmarking_1.txt
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
python train_autotalker_benchmark_models.py --dataset vizgen_merfish_mouse_liver --adata_new_name vizgen_merfish_mouse_liver_autotalker --n_neighbors_list 12 16 16 20 20 --edge_batch_size_list 256 128 128 128 128 --node_batch_size_list 32 16 16 16 16 --seeds 5 6 7 8 9 --run_index 5 6 7 8 9 10
