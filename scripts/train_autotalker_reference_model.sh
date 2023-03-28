python train_autotalker_reference_model.py \
--dataset starmap_plus_mouse_cns \
--reference_batches batch1 batch2 batch3 batch4 batch5 batch6 batch7 batch8 \
batch9 batch10 batch11 batch12 batch13 batch14 batch15 batch16 batch17 batch18 \
batch19 batch20 \
--n_neighbors 4 \
--filter_genes \
--n_hvg 2000 \
--nichenet_max_n_target_genes_per_gp 20000 \
--include_mebocost_gps \
--mebocost_species mouse \
--counts_key counts \
--log_variational \
--n_cond_embed 400 \
--n_layers_encoder 1 \
--n_epochs 40 \
--n_epochs_all_gps 20 \
--lambda_edge_recon 1000. \
--lambda_gene_expr_recon 1. \
--lambda_cond_contrastive 300. \
--cond_contrastive_thresh 0.8 \
--lambda_group_lasso 0. \
--lambda_l1_masked 0. \
--edge_batch_size 128 \
--node_batch_size 16
