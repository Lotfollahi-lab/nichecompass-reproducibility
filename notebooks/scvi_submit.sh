#!/bin/bash
#SBATCH -J SCVI_SampleBatch
#SBATCH -o logs/SCVI_Sample.out
#SBATCH -e logs/SCVI_Sample.err
#SBATCH -p gpu_p
#SBATCH -t 2-00:00:00
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --nice=10000

source ~/.bash_profile
PATH='/lustre/groups/imm01/workspace/irene.bonafonte/2023May_nichecompass/nichecompass-reproducibility'
cd $PATH
conda activate scvi-env
python notebooks/run_scvi.py -cvCat 'sample'
