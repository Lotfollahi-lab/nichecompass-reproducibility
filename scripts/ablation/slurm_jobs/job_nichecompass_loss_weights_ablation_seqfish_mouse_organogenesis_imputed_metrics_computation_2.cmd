#!/bin/bash
#SBATCH -J nichecompass_loss_weights_ablation_seqfish_mouse_organogenesis_imputed_metrics_computation_2
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_nichecompass_loss_weights_ablation_seqfish_mouse_organogenesis_imputed_metrics_computation_2.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_nichecompass_loss_weights_ablation_seqfish_mouse_organogenesis_imputed_metrics_computation_2.txt
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
python ../compute_metrics.py --task one-hop-norm_reference --file_name mlflow_summary_loss_weights_ablation_seqfish_mouse_organogenesis_imputed_2.csv --datasets seqfish_mouse_organogenesis_imputed --cell_type_keys celltype_mapped_refined --batch_keys batch --metrics gcs mlami cas clisis nasw cnmi cari casw clisi basw bgc bilisi
