#!/bin/bash
#SBATCH -J visium_mouse_brain_nichecompass_gatv2conv_sample_integration_method_benchmarking_18
#SBATCH -o ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/out_visium_mouse_brain_nichecompass_gatv2conv_sample_integration_method_benchmarking_18.txt
#SBATCH -e ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/err_visium_mouse_brain_nichecompass_gatv2conv_sample_integration_method_benchmarking_18.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/sample_integration_method_benchmarking
python ../train_nichecompass_benchmarking_models.py --adata_new_name None --n_neighbors_list 6 --edge_batch_size_list 256 --node_batch_size_list None --seeds 1 --run_index 2 --cell_type_key cell_type --filter_genes --n_svg 5000 --nichenet_keep_target_genes_ratio 1.0 --nichenet_max_n_target_genes_per_gp 250 --include_mebocost_gps --species mouse --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --dataset visium_mouse_brain --reference_batches batch1 batch2 --counts_key counts --cat_covariates_keys batch data --cat_covariates_no_edges True False --spatial_key spatial --adj_key spatial_connectivities --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label gatv2conv_sample_integration_method_benchmarking --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --n_addon_gp 10 --active_gp_thresh_ratio 0. --gene_expr_recon_dist nb --cat_covariates_embeds_injection gene_expr_decoder --cat_covariates_embeds_nums 2 2 --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gatv2conv --n_epochs 400 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 5000000. --lambda_gene_expr_recon 3000. --lambda_cat_covariates_contrastive 500000.0 --contrastive_logits_pos_ratio 0.015625 --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --lambda_l1_masked 0. --lambda_l1_addon 0. --n_sampled_neighbors 6 --timestamp_suffix _18
