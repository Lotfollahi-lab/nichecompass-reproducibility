#!/bin/sh
#SBATCH -J gpDiffTest
#SBATCH -o /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/logs/gpDiffTest.out
#SBATCH -e /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/logs/gpDiffTest.err
#SBATCH -p cpu_p
#SBATCH -t 23:00:00
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH --qos=cpu_normal
#SBATCH --nice=10000
source ~/.bash_profile
conda activate nichecompass-reproducibility
python differential_gp_testing_subgroups.py 
