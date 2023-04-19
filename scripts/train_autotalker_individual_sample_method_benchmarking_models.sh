datasets=("seqfish_mouse_organogenesis_sample2" "starmap_plus_mouse_cns_batch1" "starmap_plus_mouse_cns_subsample_20pct_batch1" "vizgen_merfish_mouse_liver" "vizgen_merfish_mouse_liver_subsample_10pct")
dataset_reference_batches=("None" "None" "None" "None" "None")
cell_type_keys=("celltype_mapped_refined" "Main_molecular_cell_type" "Main_molecular_cell_type" "Cell_Type" "Cell_Type")

len=${#datasets[@]}

for ((i=0; i<$len; i++))
do
    dataset=${datasets[i]}
    cell_type_key=${cell_type_keys[i]}
    reference_batches=${dataset_reference_batches[i]}
    
    python train_autotalker_benchmarking_models.py \
    --adata_new_name None \
    --n_neighbors_list 4 4 8 8 12 12 16 16 20 20 \
    --edge_batch_size_list 256 256 256 256 128 128 64 64 64 64 \
    --node_batch_size_list 32 32 32 32 16 16 8 8 8 8 \
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
    --model_label sample_integration_method_benchmarking \
    --active_gp_names_key autotalker_active_gp_names \
    --latent_key autotalker_latent \
    --active_gp_thresh_ratio 0.05 \
    --gene_expr_recon_dist nb \
    --cond_embed_injection gene_expr_decoder \
    --log_variational \
    --n_layers_encoder 1 \
    --conv_layer_encoder gcnconv \
    --n_epochs 40 \
    --n_epochs_all_gps 20 \
    --lr 0.001 \
    --lambda_edge_recon 1000. \
    --lambda_gene_expr_recon 1. \
    --lambda_cond_contrastive 0. \
    --contrastive_logits_ratio 0.015625 \
    --lambda_group_lasso 0. \
    --lambda_l1_masked 0.
done