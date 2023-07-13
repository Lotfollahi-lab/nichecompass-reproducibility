#!/bin/bash
#SBATCH -J nichecompass_neighbor_sampling_ablation_metrics_computation_1
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_nichecompass_neighbor_sampling_ablation_metrics_computation_1.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_nichecompass_neighbor_sampling_ablation_metrics_computation_1.txt
#SBATCH -t 12:00:00
#SBATCH -p interactive_gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive_gpu
#SBATCH --mem=64GB
#SBATCH --nice=9999
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/ablation
python ../compute_metrics.py --task neighbor_sampling_ablation --file_name mlflow_summary_neighbor_sampling_ablation_processed.csv --datasets xenium_human_breast_cancer starmap_plus_mouse_cns vizgen_merfish_human_ovarian_cancer --cell_type_keys cell_states Main_molecular_cell_type cell_type --batch_keys None None None
