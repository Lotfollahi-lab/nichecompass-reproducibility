#!/bin/bash
#SBATCH -J merfish_mouse_brain_nichecompass_postprocessing_4
#SBATCH -o ../scripts/postprocessing/slurm_jobs/logs/out_merfish_mouse_brain_nichecompass_postprocessing_4.txt
#SBATCH -e ../scripts/postprocessing/slurm_jobs/logs/err_merfish_mouse_brain_nichecompass_postprocessing_4.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=156G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/postprocessing
python ../postprocess_nichecompass_model.py --dataset merfish_mouse_brain --model_label reference --load_timestamp 21022024_194703_56 --gp_names_key nichecompass_gp_names --compute_umap --no-compute_leiden
