python train_autotalker_reference_model.py \
--dataset starmap_plus_mouse_cns \
--reference_batches batch1 batch2 batch3 \
--n_neighbors 12 \
--filter_genes \
--n_hvg 2000 \
--nichenet_max_n_target_genes_per_gp 20000 \
--include_mebocost_gps \
--mebocost_species mouse \
--counts_key counts \
--log_variational \
--n_cond_embed 475 \
--n_layers_encoder 1 \
--n_epochs 1 \
--n_epochs_all_gps 1 \
--lambda_edge_recon 10. \
--lambda_gene_expr_recon 0.1 \
--lambda_cond_contrastive 3. \
--cond_contrastive_thresh 0.8 \
--lambda_group_lasso 0. \
--lambda_l1_masked 0. \
--edge_batch_size 128 \
--node_batch_size 16
