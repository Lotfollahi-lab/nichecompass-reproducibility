#!/bin/bash
#SBATCH -J nichecompass_loss_weights_ablation_vizgen_merfish_human_ovarian_cancer_metrics_computation_3
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_nichecompass_loss_weights_ablation_vizgen_merfish_human_ovarian_cancer_metrics_computation_3.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_nichecompass_loss_weights_ablation_vizgen_merfish_human_ovarian_cancer_metrics_computation_3.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=128GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/ablation
python ../compute_metrics.py --task loss_weights_ablation --file_name mlflow_summary_loss_weights_ablation_vizgen_merfish_human_ovarian_cancer_12.csv --datasets vizgen_merfish_human_ovarian_cancer --cell_type_keys cell_type --batch_keys None --metrics gcs mlami cas clisis nasw cnmi cari casw clisi
