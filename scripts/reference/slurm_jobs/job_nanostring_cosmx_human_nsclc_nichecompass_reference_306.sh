#!/bin/bash
#SBATCH -J 306_nanostring_cosmx_human_nsclc_nichecompass_reference
#SBATCH -o logs/out_nanostring_cosmx_human_nsclc_nichecompass_reference_306.txt
#SBATCH -e logs/err_nanostring_cosmx_human_nsclc_nichecompass_reference_306.txt
#SBATCH -t 24:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_normal
#SBATCH --mem=128G
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass-reproducibility
cd /lustre/groups/imm01/workspace/irene.bonafonte/Projects/2023May_nichecompass/nichecompass-reproducibility/scripts/reference

python ../train_nichecompass_reference_model.py --dataset nanostring_cosmx_human_nsclc --timestamp_suffix _306 --reference_batches batch1 batch2 batch4 batch5 batch6 batch7  --cat_covariates_embeds_nums 3 30 5 --contrastive_logits_pos_ratio 0.03125 --lambda_l1_addon 1000 --lambda_l1_masked 0 --active_gp_thresh_ratio 0.03 --lambda_cat_covariates_contrastive 1000000.0 --n_neighbors 4 --no-filter_genes --nichenet_keep_target_genes_ratio 1.0 --nichenet_max_n_target_genes_per_gp 250 --include_mebocost_gps --species human --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --counts_key counts --cat_covariates_keys batch fov patient --cat_covariates_no_edges True False True --spatial_key spatial --adj_key spatial_connectivities --mapping_entity_key mapping_entity --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label reference --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --gene_expr_recon_dist nb --cat_covariates_embeds_injection gene_expr_decoder --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --conv_layer_encoder gatv2conv --n_epochs 400 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 5000000. --lambda_gene_expr_recon 3000. --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --n_sampled_neighbors 4  --edge_batch_size 512 --node_batch_size 256 --n_hidden_encoder 960 --n_addon_gp 100