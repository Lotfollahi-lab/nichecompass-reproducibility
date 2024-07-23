#!/bin/bash
#SBATCH -J graphst_single_sample_method_benchmarking_sim1_1105genes_10000locs_strongincrements_metrics_computation_1
#SBATCH -o ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/out_graphst_single_sample_method_benchmarking_sim1_1105genes_10000locs_strongincrements_metrics_computation_1.txt
#SBATCH -e ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/err_graphst_single_sample_method_benchmarking_sim1_1105genes_10000locs_strongincrements_metrics_computation_1.txt
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
python ../compute_benchmarking_metrics.py --dataset sim1_1105genes_10000locs_strongincrements --task single_sample_method_benchmarking --file_name sim1_1105genes_10000locs_strongincrements_graphst.h5ad --cell_type_key cell_type --niche_type_key niche_type --batch_key None --latent_key graphst_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi --include_sdmbench
