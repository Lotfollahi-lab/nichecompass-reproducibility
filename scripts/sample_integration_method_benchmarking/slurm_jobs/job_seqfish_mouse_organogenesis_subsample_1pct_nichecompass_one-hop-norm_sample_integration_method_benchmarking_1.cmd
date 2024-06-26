#!/bin/bash
#SBATCH -J seqfish_mouse_organogenesis_subsample_1pct_nichecompass_one-hop-norm_sample_integration_method_benchmarking_1
#SBATCH -o ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/out_seqfish_mouse_organogenesis_subsample_1pct_nichecompass_one-hop-norm_sample_integration_method_benchmarking_1.txt
#SBATCH -e ../scripts/sample_integration_method_benchmarking/slurm_jobs/logs/err_seqfish_mouse_organogenesis_subsample_1pct_nichecompass_one-hop-norm_sample_integration_method_benchmarking_1.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/sample_integration_method_benchmarking
python ../train_nichecompass_benchmarking_models.py --adata_new_name None --n_neighbors_list 4 4 8 8 12 12 16 16 20 20 --edge_batch_size_list 16384 16384 16384 16384 16384 16384 16384 16384 16384 16384 --node_batch_size_list None None None None None None None None None None --seeds 0 1 2 3 4 5 6 7 8 9 --run_index 1 2 3 4 5 6 7 8 9 10 --cell_type_key celltype_mapped_refined --no-filter_genes --nichenet_keep_target_genes_ratio 1. --nichenet_max_n_target_genes_per_gp 250 --include_mebocost_gps --species mouse --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --dataset seqfish_mouse_organogenesis_subsample_1pct --reference_batches batch1 batch2 batch3 batch4 batch5 batch6 --counts_key counts --cat_covariates_keys batch --spatial_key spatial --adj_key spatial_connectivities --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label one-hop-norm_sample_integration_method_benchmarking --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --active_gp_thresh_ratio 0.05 --gene_expr_recon_dist nb --cat_covariates_embeds_injection gene_expr_decoder --cat_covariates_embeds_nums 6 --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gcnconv --n_epochs 100 --n_epochs_all_gps 25 --n_epochs_no_cat_covariates_contrastive 0 --lr 0.001 --lambda_edge_recon 500000. --lambda_gene_expr_recon 300. --lambda_cat_covariates_contrastive 0.0 --contrastive_logits_pos_ratio 0.0 --contrastive_logits_neg_ratio 0.0 --lambda_group_lasso 0. --lambda_l1_masked 0. --timestamp_suffix _1
