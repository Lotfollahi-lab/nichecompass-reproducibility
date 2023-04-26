python train_autotalker_reference_model.py \
--dataset seqfish_mouse_organogenesis_imputed \
--reference_batches batch1 batch2 batch3 batch4 batch5 batch6 \
--n_neighbors 12 \
--filter_genes \
--n_hvg 2000 \
--nichenet_keep_target_genes_ratio 0.01 \
--nichenet_max_n_target_genes_per_gp 25344 \
--include_mebocost_gps \
--mebocost_species mouse \
--gp_filter_mode subset \
--combine_overlap_gps \
--overlap_thresh_source_genes 0.9 \
--overlap_thresh_target_genes 0.9 \
--overlap_thresh_genes 0.9 \
--counts_key log_normalized_counts \
--condition_key batch \
--spatial_key spatial \
--adj_key spatial_connectivities \
--mapping_entity_key mapping_entity \
--gp_targets_mask_key autotalker_gp_targets \
--gp_sources_mask_key autotalker_gp_sources \
--gp_names_key autotalker_gp_names \
--model_label reference \
--active_gp_names_key autotalker_active_gp_names \
--latent_key autotalker_latent \
--active_gp_thresh_ratio 0.05 \
--gene_expr_recon_dist nb \
--cond_embed_injection gene_expr_decoder \
--no-log_variational \
--node_label_method one-hop-norm \
--n_layers_encoder 1 \
--conv_layer_encoder gcnconv \
--n_epochs 100 \
--n_epochs_all_gps 25 \
--n_epochs_no_cond_contrastive 0 \
--lambda_edge_recon 500000. \
--lambda_gene_expr_recon 100. \
--lambda_cond_contrastive 0. \
--contrastive_logits_ratio 0.125 \
--lambda_group_lasso 0. \
--lambda_l1_masked 10. \
--edge_batch_size 16384 \
--node_batch_size None