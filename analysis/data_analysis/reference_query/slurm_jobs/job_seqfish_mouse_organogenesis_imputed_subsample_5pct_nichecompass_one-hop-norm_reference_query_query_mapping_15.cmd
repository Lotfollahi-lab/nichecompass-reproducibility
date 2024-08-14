#!/bin/bash
#SBATCH -J seqfish_mouse_organogenesis_imputed_subsample_5pct_nichecompass_one-hop-norm_reference_query_query_mapping_15
#SBATCH -o ./data_analysis/reference_query/slurm_jobs/logs/out_seqfish_mouse_organogenesis_imputed_subsample_5pct_nichecompass_one-hop-norm_reference_query_query_mapping_15.txt
#SBATCH -e ./data_analysis/reference_query/slurm_jobs/logs/err_seqfish_mouse_organogenesis_imputed_subsample_5pct_nichecompass_one-hop-norm_reference_query_query_mapping_15.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/analysis/data_analysis/reference_query
python ../map_query_on_nichecompass_reference_model.py --reference_dataset seqfish_mouse_organogenesis_imputed_subsample_5pct --query_dataset seqfish_mouse_organogenesis_imputed --query_batches batch5 batch6 --n_neighbors 8 --spatial_key spatial --mapping_entity_key mapping_entity --gp_names_key nichecompass_gp_names --reference_model_label reference --load_timestamp 08082024_140042_15 --query_model_label query --reference_query_model_label reference_query --n_epochs 100 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 5000000. --lambda_gene_expr_recon 3000. --lambda_cat_covariates_contrastive 0.0 --contrastive_logits_pos_ratio 0.0 --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --lambda_l1_masked 0. --edge_batch_size 256 --node_batch_size None
