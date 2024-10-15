#!/bin/bash
#SBATCH -J starmap_plus_mouse_cns_nichecompass_loss_extended_ablation_42
#SBATCH -o ./ablation/slurm_jobs/logs/out_starmap_plus_mouse_cns_nichecompass_loss_extended_ablation_42.txt
#SBATCH -e ./ablation/slurm_jobs/logs/err_starmap_plus_mouse_cns_nichecompass_loss_extended_ablation_42.txt
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
cd /home/aih/sebastian.birk_new/workspace/projects/nichecompass-reproducibility/analysis/ablation
python .././data_analysis/train_nichecompass_reference_model.py --dataset starmap_plus_mouse_cns --reference_batches batch1 --n_neighbors 6 --no-filter_genes --nichenet_keep_target_genes_ratio 1. --nichenet_max_n_target_genes_per_gp 250 --include_mebocost_gps --species mouse --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --counts_key counts --spatial_key spatial --adj_key spatial_connectivities --mapping_entity_key mapping_entity --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label loss_extended_ablation --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --active_gp_thresh_ratio 0.01 --gene_expr_recon_dist nb --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gatv2conv --n_epochs 400 --n_epochs_all_gps 25 --lr 0.001 --lambda_edge_recon 500000.0 --lambda_gene_expr_recon 300.0 --lambda_group_lasso 0. --lambda_l1_masked 30.0 --lambda_l1_addon 30.0 --l1_targets_categories target_gene --l1_sources_categories None --edge_batch_size 256 --node_batch_size None --n_sampled_neighbors 4 --seed 1 --use_new_gp_mask --timestamp_suffix _42
