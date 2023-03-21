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
--no-log_variational \
--n_epochs 40 \
--n_epochs_all_gps 20 \
--lambda_group_lasso 0. \
--lambda_l1_masked 0. \
--edge_batch_size 256 \
--node_batch_size 32
