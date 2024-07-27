#!/bin/bash
#SBATCH -J starmap_mouse_mpfc_nichecompass_gcnconv_single_sample_method_benchmarking_1
#SBATCH -o ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/out_starmap_mouse_mpfc_nichecompass_gcnconv_single_sample_method_benchmarking_1.txt
#SBATCH -e ./benchmarking/single_sample_method_benchmarking/slurm_jobs/logs/err_starmap_mouse_mpfc_nichecompass_gcnconv_single_sample_method_benchmarking_1.txt
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
cd /home/aih/sebastian.birk/workspace/projects/nichecompass-reproducibility/analysis/benchmarking/single_sample_method_benchmarking
python ../train_nichecompass_benchmarking_models.py --adata_new_name None --n_neighbors_list 16 16 12 12 8 8 4 4 --edge_batch_size_list 131072 131072 131072 131072 131072 131072 131072 131072 --node_batch_size_list None None None None None None None None --seeds 7 6 5 4 3 2 1 0 --run_index 8 7 6 5 4 3 2 1 --cell_type_key cell_type --nichenet_keep_target_genes_ratio 1.0 --nichenet_max_n_target_genes_per_gp 250 --include_mebocost_gps --species mouse --gp_filter_mode subset --combine_overlap_gps --overlap_thresh_source_genes 0.9 --overlap_thresh_target_genes 0.9 --overlap_thresh_genes 0.9 --dataset starmap_mouse_mpfc --reference_batches None --counts_key counts --spatial_key spatial --adj_key spatial_connectivities --mapping_entity_key mapping_entity --no-filter_genes --gp_targets_mask_key nichecompass_gp_targets --gp_sources_mask_key nichecompass_gp_sources --gp_names_key nichecompass_gp_names --model_label gcnconv_single_sample_method_benchmarking --active_gp_names_key nichecompass_active_gp_names --latent_key nichecompass_latent --n_addon_gp 10 --active_gp_thresh_ratio 0. --gene_expr_recon_dist nb --log_variational --node_label_method one-hop-norm --n_layers_encoder 1 --n_hidden_encoder None --conv_layer_encoder gcnconv --n_epochs 400 --n_epochs_all_gps 25 --lr 0.001 --lambda_edge_recon 5000000. --lambda_gene_expr_recon 3000. --lambda_group_lasso 0. --lambda_l1_masked 0. --lambda_l1_addon 0. --n_sampled_neighbors 4 --use_new_gp_mask --timestamp_suffix _1
