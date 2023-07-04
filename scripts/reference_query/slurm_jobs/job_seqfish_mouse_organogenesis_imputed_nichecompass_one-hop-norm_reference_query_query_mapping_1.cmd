#!/bin/bash
#SBATCH -J seqfish_mouse_organogenesis_imputed_nichecompass_one-hop-norm_reference_query_query_mapping_1
#SBATCH -o ../scripts/reference_query/slurm_jobs/logs/out_seqfish_mouse_organogenesis_imputed_nichecompass_one-hop-norm_reference_query_query_mapping_1.txt
#SBATCH -e ../scripts/reference_query/slurm_jobs/logs/err_seqfish_mouse_organogenesis_imputed_nichecompass_one-hop-norm_reference_query_query_mapping_1.txt
#SBATCH -t 12:00:00
#SBATCH -p interactive_gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive_gpu
#SBATCH --mem=64GB
#SBATCH --nice=9999
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference_query
python ../map_query_on_nichecompass_reference_model.py --dataset seqfish_mouse_organogenesis_imputed --query_batches batch5 batch6 --n_neighbors 12 --spatial_key spatial --mapping_entity_key mapping_entity --gp_names_key nichecompass_gp_names --reference_model_label one-hop-norm_reference_query_reference_only --load_timestamp 01072023_165203_1 --query_model_label one-hop-norm_reference_query_query_only --reference_query_model_label one-hop-norm_reference_query_query_mapping --n_epochs 100 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 500000. --lambda_gene_expr_recon 300. --lambda_cat_covariates_contrastive 0.0 --contrastive_logits_pos_ratio 0.0 --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --lambda_l1_masked 5.0 --edge_batch_size 4096 --node_batch_size None
