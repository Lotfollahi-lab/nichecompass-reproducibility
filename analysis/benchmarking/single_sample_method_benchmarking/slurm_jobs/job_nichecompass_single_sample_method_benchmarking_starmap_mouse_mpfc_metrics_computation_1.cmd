#!/bin/bash
#SBATCH -J nichecompass_single_sample_method_benchmarking_starmap_mouse_mpfc_metrics_computation_1
#SBATCH -o ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/out_nichecompass_single_sample_method_benchmarking_starmap_mouse_mpfc_metrics_computation_1.txt
#SBATCH -e ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/err_nichecompass_single_sample_method_benchmarking_starmap_mouse_mpfc_metrics_computation_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=156G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/analysis/benchmarking/single_sample_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset starmap_mouse_mpfc --task single_sample_method_benchmarking --file_name starmap_mouse_mpfc_nichecompass_gatv2conv.h5ad --cell_type_key cell_type --niche_type_key niche_type --batch_key None --latent_key nichecompass_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi --include_sdmbench
