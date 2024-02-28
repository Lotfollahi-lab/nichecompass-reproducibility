#!/bin/bash
#SBATCH -J merfish_mouse_brain_nichecompass_postprocessing_7
#SBATCH -o ../scripts/postprocessing/slurm_jobs/logs/out_merfish_mouse_brain_nichecompass_postprocessing_7.txt
#SBATCH -e ../scripts/postprocessing/slurm_jobs/logs/err_merfish_mouse_brain_nichecompass_postprocessing_7.txt
#SBATCH -t 48:00:00
#SBATCH -p cpu_p
#SBATCH -c 6
#SBATCH --qos=cpu_normal
#SBATCH --mem=300G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/postprocessing
python ../postprocess_nichecompass_model.py --dataset merfish_mouse_brain --model_label reference --load_timestamp 24022024_170907_82 --gp_names_key nichecompass_gp_names --compute_latent --compute_pca --compute_knn_graph --compute_umap --compute_leiden
