#!/bin/bash
#SBATCH -J banksy_single_sample_method_benchmarking_vizgen_merfish_mouse_liver_subsample_1pct_metrics_computation_1
#SBATCH -o ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/out_banksy_single_sample_method_benchmarking_vizgen_merfish_mouse_liver_subsample_1pct_metrics_computation_1.txt
#SBATCH -e ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/err_banksy_single_sample_method_benchmarking_vizgen_merfish_mouse_liver_subsample_1pct_metrics_computation_1.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/analysis/benchmarking/single_sample_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset vizgen_merfish_mouse_liver_subsample_1pct --task single_sample_method_benchmarking --file_name vizgen_merfish_mouse_liver_subsample_1pct_banksy.h5ad --cell_type_key cell_type --batch_key None --latent_key banksy_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi
