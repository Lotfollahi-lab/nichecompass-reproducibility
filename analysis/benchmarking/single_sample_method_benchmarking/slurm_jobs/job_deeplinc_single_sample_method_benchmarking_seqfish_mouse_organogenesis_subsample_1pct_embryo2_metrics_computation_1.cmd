#!/bin/bash
#SBATCH -J deeplinc_single_sample_method_benchmarking_seqfish_mouse_organogenesis_subsample_1pct_embryo2_metrics_computation_1
#SBATCH -o ../scripts/single_sample_method_benchmarking/slurm_jobs/logs/out_deeplinc_single_sample_method_benchmarking_seqfish_mouse_organogenesis_subsample_1pct_embryo2_metrics_computation_1.txt
#SBATCH -e ../scripts/single_sample_method_benchmarking/slurm_jobs/logs/err_deeplinc_single_sample_method_benchmarking_seqfish_mouse_organogenesis_subsample_1pct_embryo2_metrics_computation_1.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset seqfish_mouse_organogenesis_subsample_1pct_embryo2 --task single_sample_method_benchmarking --file_name seqfish_mouse_organogenesis_subsample_1pct_embryo2_deeplinc.h5ad --cell_type_key cell_type --batch_key None --latent_key deeplinc_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi
