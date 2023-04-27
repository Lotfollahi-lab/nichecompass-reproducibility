datasets=("starmap_plus_mouse_cns")
dataset_reference_batches=("batch1 batch2 batch3 batch4 batch5 batch6 batch7 batch8 batch9 batch10 batch11 batch12 batch13 batch14 batch15 batch16 batch17 batch18 batch19 batch20")
cell_type_keys=("Main_molecular_cell_type")

len=${#datasets[@]}

python waiting_script.py &&

for ((i=0; i<$len; i++))
do
    dataset=${datasets[i]}
    cell_type_key=${cell_type_keys[i]}
    reference_batches=${dataset_reference_batches[i]}
    
    python train_autotalker_benchmarking_models.py \
    --adata_new_name None \
    --n_neighbors_list 4 4 8 8 12 12 16 16 20 20 \
    --edge_batch_size_list 2048 2048 1024 1024 512 512 256 256 128 128 \
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
    --gp_targets_mask_key autotalker_gp_targets \
    --gp_sources_mask_key autotalker_gp_sources \
    --gp_names_key autotalker_gp_names \
    --model_label one-hop-attention_2_sample_integration_method_benchmarking \
    --active_gp_names_key autotalker_active_gp_names \
    --latent_key autotalker_latent \
    --active_gp_thresh_ratio 0.05 \
    --gene_expr_recon_dist nb \
    --cond_embed_injection gene_expr_decoder \
    --log_variational \
    --node_label_method one-hop-attention \
    --n_layers_encoder 1 \
    --conv_layer_encoder gcnconv \
    --n_epochs 100 \
    --n_epochs_all_gps 25 \
    --n_epochs_no_cond_contrastive 5 \
    --lr 0.001 \
    --lambda_edge_recon 500000. \
    --lambda_gene_expr_recon 100. \
    --lambda_cond_contrastive 100000. \
    --contrastive_logits_ratio 0.015625 \
    --lambda_group_lasso 0. \
    --lambda_l1_masked 0.
done