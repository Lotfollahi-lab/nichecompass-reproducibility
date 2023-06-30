#!/bin/bash
#SBATCH -J nichecompass_loss_weights_ablation_metrics_computation_1
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_nichecompass_loss_weights_ablation_metrics_computation_1.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_nichecompass_loss_weights_ablation_metrics_computation_1.txt
#SBATCH -t 12:00:00
#SBATCH -p interactive_gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive_gpu
#SBATCH --mem=64GB
#SBATCH --nice=9999
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/ablation
python ../compute_metrics.py --ablation_task loss_weights --datasets xenium_human_breast_cancer starmap_plus_mouse_cns --cell_type_keys cell_states Main_molecular_cell_type --condition_keys None None --experiment_ids 3 4
