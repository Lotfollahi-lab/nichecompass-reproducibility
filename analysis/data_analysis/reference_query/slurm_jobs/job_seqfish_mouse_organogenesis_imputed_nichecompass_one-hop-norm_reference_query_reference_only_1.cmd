#!/bin/bash
#SBATCH -J seqfish_mouse_organogenesis_imputed_nichecompass_one-hop-norm_reference_query_reference_only_1
#SBATCH -o ../scripts/reference_query/slurm_jobs/logs/out_seqfish_mouse_organogenesis_imputed_nichecompass_one-hop-norm_reference_query_reference_only_1.txt
#SBATCH -e ../scripts/reference_query/slurm_jobs/logs/err_seqfish_mouse_organogenesis_imputed_nichecompass_one-hop-norm_reference_query_reference_only_1.txt
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
python ../train_nichecompass_reference_model.py --dataset seqfish_mouse_organogenesis_imputed --reference_batches batch1 batch2 batch3 batch4 --n_neighbors 12 --filter_genes --n_hvg 3000 --nichenet_keep_target_genes_ratio 0.1 --nichenet_max_n_target_genes_per_gp 250 --include_mebocost_gps --species mouse --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --counts_key counts --condition_key batch --cat_covariates_keys batch --cat_covariates_no_edges True --spatial_key spatial --adj_key spatial_connectivities --mapping_entity_key mapping_entity --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label one-hop-norm_reference_query_reference_only --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --active_gp_thresh_ratio 0.05 --gene_expr_recon_dist nb --cat_covariates_embeds_injection gene_expr_decoder --cat_covariates_embeds_nums 3 --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gcnconv --n_epochs 100 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 500000. --lambda_gene_expr_recon 300. --lambda_cat_covariates_contrastive 0.0 --contrastive_logits_pos_ratio 0.0 --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --lambda_l1_masked 5.0 --edge_batch_size 4096 --node_batch_size None --timestamp_suffix _1
