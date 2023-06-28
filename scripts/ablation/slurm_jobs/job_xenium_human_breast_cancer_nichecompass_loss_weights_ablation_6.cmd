#!/bin/bash
#SBATCH -J xenium_human_breast_cancer_nichecompass_loss_weights_ablation_6
#SBATCH -o ../scripts/ablation/slurm_jobs/logs/out_xenium_human_breast_cancer_nichecompass_loss_weights_ablation_6.txt
#SBATCH -e ../scripts/ablation/slurm_jobs/logs/err_xenium_human_breast_cancer_nichecompass_loss_weights_ablation_6.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/scripts/ablation
python ../train_nichecompass_reference_model.py --dataset xenium_human_breast_cancer --reference_batches batch1 --n_neighbors 4 --no-filter_genes --nichenet_keep_target_genes_ratio 1. --nichenet_max_n_target_genes_per_gp 250 --add_fc_gps_instead_of_gp_dict_gps --include_mebocost_gps --species human --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --counts_key counts --spatial_key spatial --adj_key spatial_connectivities --mapping_entity_key mapping_entity --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label loss_weights_ablation --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --active_gp_thresh_ratio 0. --gene_expr_recon_dist nb --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gcnconv --n_epochs 100 --n_epochs_all_gps 25 --lr 0.001 --lambda_edge_recon 0 --lambda_gene_expr_recon 30 --lambda_group_lasso 0. --lambda_l1_masked 0. --edge_batch_size 512 --node_batch_size None --timestamp_suffix _6
