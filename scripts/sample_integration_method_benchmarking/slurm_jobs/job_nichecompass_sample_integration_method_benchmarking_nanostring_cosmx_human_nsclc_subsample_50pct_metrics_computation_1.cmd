#!/bin/bash
#SBATCH -J nichecompass_sample_integration_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_50pct_metrics_computation_1
#SBATCH -o ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/out_nichecompass_sample_integration_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_50pct_metrics_computation_1.txt
#SBATCH -e ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/err_nichecompass_sample_integration_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_50pct_metrics_computation_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=128G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/sample_integration_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset nanostring_cosmx_human_nsclc_subsample_50pct --task sample_integration_method_benchmarking --file_name nanostring_cosmx_human_nsclc_subsample_50pct_nichecompass_gatv2conv_fov.h5ad --cell_type_key cell_type --batch_key batch --latent_key nichecompass_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi basw bgc blisi
