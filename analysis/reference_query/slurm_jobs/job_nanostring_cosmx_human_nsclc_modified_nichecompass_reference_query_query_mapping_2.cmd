#!/bin/bash
#SBATCH -J nanostring_cosmx_human_nsclc_modified_nichecompass_reference_query_query_mapping_2
#SBATCH -o ../scripts/reference_query/slurm_jobs/logs/out_nanostring_cosmx_human_nsclc_modified_nichecompass_reference_query_query_mapping_2.txt
#SBATCH -e ../scripts/reference_query/slurm_jobs/logs/err_nanostring_cosmx_human_nsclc_modified_nichecompass_reference_query_query_mapping_2.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=128G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-test
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference_query
python ../map_query_on_nichecompass_reference_model.py --dataset nanostring_cosmx_human_nsclc_modified --query_batches batch3 --n_neighbors 4 --spatial_key spatial --mapping_entity_key mapping_entity --gp_names_key nichecompass_gp_names --reference_model_label reference --load_timestamp 01092023_182316_2 --query_model_label query --reference_query_model_label reference_query_mapping --n_epochs 400 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 5000000. --lambda_gene_expr_recon 3000. --lambda_cat_covariates_contrastive 0.0 --contrastive_logits_pos_ratio 0.0 --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --lambda_l1_masked 0.0 --edge_batch_size 512 --node_batch_size None --n_sampled_neighbors 4
