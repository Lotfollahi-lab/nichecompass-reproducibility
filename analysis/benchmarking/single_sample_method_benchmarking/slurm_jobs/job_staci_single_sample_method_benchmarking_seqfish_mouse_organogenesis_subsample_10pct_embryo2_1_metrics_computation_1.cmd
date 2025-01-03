#!/bin/bash
#SBATCH -J staci_single_sample_method_benchmarking_seqfish_mouse_organogenesis_subsample_10pct_embryo2_1_metrics_computation_1
#SBATCH -o ../scripts/single_sample_method_benchmarking/slurm_jobs/logs/out_staci_single_sample_method_benchmarking_seqfish_mouse_organogenesis_subsample_10pct_embryo2_1_metrics_computation_1.txt
#SBATCH -e ../scripts/single_sample_method_benchmarking/slurm_jobs/logs/err_staci_single_sample_method_benchmarking_seqfish_mouse_organogenesis_subsample_10pct_embryo2_1_metrics_computation_1.txt
#SBATCH -t 48:00:00
#SBATCH -p cpu_p
#SBATCH -c 6
#SBATCH --qos=cpu_normal
#SBATCH --mem=300G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/single_sample_method_benchmarking
python ../compute_benchmarking_metrics.py --dataset seqfish_mouse_organogenesis_subsample_10pct_embryo2 --task single_sample_method_benchmarking --file_name seqfish_mouse_organogenesis_subsample_10pct_embryo2_staci_1.h5ad --cell_type_key cell_type --batch_key None --latent_key staci_latent --metrics gcs mlami cas clisis nasw cnmi cari casw clisi --run_numbers 1
