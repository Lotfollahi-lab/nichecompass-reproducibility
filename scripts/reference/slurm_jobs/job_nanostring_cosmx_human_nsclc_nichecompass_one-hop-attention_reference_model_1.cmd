#!/bin/bash
#SBATCH -J nanostring_cosmx_human_nsclc_nichecompass_one-hop-attention_reference_model_1
#SBATCH -o ../scripts/reference/slurm_jobs/logs/out_nanostring_cosmx_human_nsclc_nichecompass_one-hop-attention_reference_model_1.txt
#SBATCH -e ../scripts/reference/slurm_jobs/logs/err_nanostring_cosmx_human_nsclc_nichecompass_one-hop-attention_reference_model_1.txt
#SBATCH -t 48:00:00
#SBATCH -p gpu_p
#SBATCH -c 6
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu
#SBATCH --mem=64GB
#SBATCH --nice=10000
source $HOME/.bashrc
conda activate nichecompass
cd /
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/reference
python ../train_nichecompass_reference_model.py --dataset nanostring_cosmx_human_nsclc --reference_batches batch1 batch2 batch3 batch4 batch5 batch6 batch7 batch8 --n_neighbors 12 --no-filter_genes --nichenet_keep_target_genes_ratio 0.01 --nichenet_max_n_target_genes_per_gp 1000 --include_mebocost_gps --mebocost_species human --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --counts_key counts --condition_key batch --spatial_key spatial --adj_key spatial_connectivities --mapping_entity_key mapping_entity --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label one-hop-attention_reference --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --active_gp_thresh_ratio 0.1 --gene_expr_recon_dist nb --cond_embed_injection gene_expr_decoder --n_cond_embed 8 --log_variational --node_label_method one-hop-attention --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gcnconv --n_epochs 100 --n_epochs_all_gps 25 --n_epochs_no_cond_contrastive 0 --lr 0.001 --lambda_edge_recon 500000. --lambda_gene_expr_recon 300. --lambda_cond_contrastive 0. --contrastive_logits_ratio 0. --lambda_group_lasso 0. --lambda_l1_masked 5. --edge_batch_size 1024 --node_batch_size None --timestamp_suffix _1
