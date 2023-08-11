#!/bin/bash
#SBATCH -J nichecompass_loss_weights_ablation_starmap_plus_mouse_cns_metrics_computation_2
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_nichecompass_loss_weights_ablation_starmap_plus_mouse_cns_metrics_computation_2.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_nichecompass_loss_weights_ablation_starmap_plus_mouse_cns_metrics_computation_2.txt
#SBATCH -t 48:00:00
#SBATCH -p cpu_p
#SBATCH -c 6
#SBATCH --mem=128GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/ablation
python ../compute_metrics.py --task loss_weights_ablation --file_name mlflow_summary_loss_weights_ablation_starmap_plus_mouse_cns_50.csv --datasets starmap_plus_mouse_cns --cell_type_keys Main_molecular_cell_type --batch_keys None --metrics gcs mlami cas clisis nasw cnmi cari casw clisi
