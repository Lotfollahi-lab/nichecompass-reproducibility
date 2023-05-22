#!/bin/bash
#SBATCH -J vizgen_merfish_human_ovarian_cancer_nichecompass_one-hop-attention_reference_model
#SBATCH -o /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference/slurm_jobs/logs/out_vizgen_merfish_human_ovarian_cancer_nichecompass_one-hop-attention_reference_model.txt
#SBATCH -e /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference/slurm_jobs/logs/err_vizgen_merfish_human_ovarian_cancer_nichecompass_one-hop-attention_reference_model.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference
bash train_vizgen_merfish_human_ovarian_cancer_nichecompass_one-hop-attention_reference_model.sh
