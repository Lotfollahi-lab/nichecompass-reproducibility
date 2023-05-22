datasets=("starmap_plus_mouse_cns_batch1")
dataset_reference_batches=("None")
cell_type_keys=("Main_molecular_cell_type")

len=${#datasets[@]}

for ((i=0; i<$len; i++))
do
    dataset=${datasets[i]}
    cell_type_key=${cell_type_keys[i]}
    reference_batches=${dataset_reference_batches[i]}
    
    python ../wait.py
    python ../train_nichecompass_benchmarking_models.py \
    --adata_new_name None \
    --n_neighbors_list 4 4 8 8 12 12 16 16 20 20 \
    --edge_batch_size_list 4096 4096 2048 2048 1024 1024 512 512 256 256 \
    --node_batch_size_list None None None None None None None None None None \
    --seeds 0 1 2 3 4 5 6 7 8 9 \
    --run_index 1 2 3 4 5 6 7 8 9 10 \
    --cell_type_key $cell_type_key \
    --nichenet_keep_target_genes_ratio 0.01 \
    --nichenet_max_n_target_genes_per_gp 25344 \
    --include_mebocost_gps \
    --mebocost_species mouse \
    --gp_filter_mode subset \
    --combine_overlap_gps \
    --overlap_thresh_source_genes 0.9 \
    --overlap_thresh_target_genes 0.9 \
    --overlap_thresh_genes 0.9 \
    --dataset $dataset \
    --reference_batches $reference_batches \
    --counts_key counts \
    --condition_key batch \
    --spatial_key spatial \
    --adj_key spatial_connectivities \
    --mapping_entity_key mapping_entity \
    --no-filter_genes \
    --gp_targets_mask_key nichecompass_gp_targets \
    --gp_sources_mask_key nichecompass_gp_sources \
    --gp_names_key nichecompass_gp_names \
    --model_label one-hop-attention_single_sample_method_benchmarking \
    --active_gp_names_key nichecompass_active_gp_names \
    --latent_key nichecompass_latent \
    --active_gp_thresh_ratio 0.1 \
    --gene_expr_recon_dist nb \
    --cond_embed_injection gene_expr_decoder \
    --n_cond_embed None \
    --log_variational \
    --node_label_method one-hop-attention \
    --n_layers_encoder 1 \
    --n_hidden_encoder None \
    --conv_layer_encoder gcnconv \
    --n_epochs 100 \
    --n_epochs_all_gps 25 \
    --n_epochs_no_cond_contrastive 0 \
    --lr 0.001 \
    --lambda_edge_recon 500000. \
    --lambda_gene_expr_recon 300. \
    --lambda_cond_contrastive 0. \
    --contrastive_logits_ratio 0. \
    --lambda_group_lasso 0. \
    --lambda_l1_masked 0.
done