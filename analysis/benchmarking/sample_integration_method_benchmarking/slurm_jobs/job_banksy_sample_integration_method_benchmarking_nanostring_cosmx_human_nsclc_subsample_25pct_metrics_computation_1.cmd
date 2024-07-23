#!/bin/bash
#SBATCH -J banksy_sample_integration_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_25pct_metrics_computation_1
#SBATCH -o ./benchmarking/sample_integration_method_benchmarking/slurm_jobs/logs/out_banksy_sample_integration_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_25pct_metrics_computation_1.txt
#SBATCH -e ./benchmarking/sample_integration_method_benchmarking/slurm_jobs/logs/err_banksy_sample_integration_method_benchmarking_nanostring_cosmx_human_nsclc_subsample_25pct_metrics_computation_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=156G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/analysis/benchmarking/sample_integration_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset nanostring_cosmx_human_nsclc_subsample_25pct --task sample_integration_method_benchmarking --file_name nanostring_cosmx_human_nsclc_subsample_25pct_banksy.h5ad --cell_type_key cell_type --batch_key batch --batches batch1 batch2 batch3 --latent_key banksy_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi basw blisi kbet pcr
