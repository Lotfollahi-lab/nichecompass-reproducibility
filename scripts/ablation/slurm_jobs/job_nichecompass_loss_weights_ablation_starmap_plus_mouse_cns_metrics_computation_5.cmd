#!/bin/bash
#SBATCH -J nichecompass_loss_weights_ablation_starmap_plus_mouse_cns_metrics_computation_5
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_nichecompass_loss_weights_ablation_starmap_plus_mouse_cns_metrics_computation_5.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_nichecompass_loss_weights_ablation_starmap_plus_mouse_cns_metrics_computation_5.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/ablation
python ../compute_metrics.py --task loss_weights_ablation --file_name mlflow_summary_loss_weights_ablation_starmap_plus_mouse_cns_100.csv --datasets starmap_plus_mouse_cns --cell_type_keys Main_molecular_cell_type --batch_keys None --metrics gcs nasw cnmi
