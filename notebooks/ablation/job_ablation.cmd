#!/bin/bash
#SBATCH -J ablation_metrics
#SBATCH -o ../notebooks/ablation/out_ablation.txt
#SBATCH -e ../notebooks/ablation/err_ablation.txt
#SBATCH -t 48:00:00
#SBATCH -p cpu_p
#SBATCH -c 6
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/notebooks/ablation
python ablation.py
