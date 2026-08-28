# Shared configuration for the Xenium human breast cancer run with human PPI
# gene programs, pinned to xenium_human_breast_cancer_humanppi_gps.ipynb.
#
# Sourced by both the single device warm-up job and the multi-GPU job, so that
# the two cannot drift apart. Nothing here is specific to one device count or
# to one scheduler: ´--multi_gpu´, ´--n_epochs´ and ´--model_label´ are added by
# the caller.
#
# Data is read from datasets/st_data/gold/{dataset}_{batch}.h5ad, the same
# location the notebook reads from.

NICHECOMPASS_ARGS=(
    --dataset xenium_human_breast_cancer
    --reference_batches batch1 batch2
    --counts_key counts
    --spatial_key spatial
    --n_neighbors 8
    --species human

    # Prior gene programs: the human PPI resource only, as in the notebook
    --no-include_omnipath_gps
    --no-include_nichenet_gps
    --no-include_mebocost_gps
    --no-include_collectri_gps
    --include_humanppi_gps
    --humanppi_precision 80
    --humanppi_program_type intercellular
    --humanppi_ambiguous_locality extracellular
    --humanppi_unresolved_locality exclude
    --humanppi_use_topology
    --humanppi_detect_cis_complexes
    --humanppi_min_extracellular_domain_length 30
    --humanppi_orient_juxtacrine_gps
    --no-humanppi_symmetric_juxtacrine_gps
    --humanppi_filter_ig_tcr_segments
    --humanppi_filter_paralog_cross_pairs
    --humanppi_min_rf_prob None
    --humanppi_min_af_prob None

    # No combining, and keep only programs that still encode an interaction
    --gp_filter_mode none
    --min_genes_per_gp 2
    --min_source_genes_per_gp 1
    --min_target_genes_per_gp 1

    # Model architecture
    --conv_layer_encoder gatv2conv
    --active_gp_thresh_ratio 0.01
    --n_addon_gp 100

    # Batch as a categorical covariate, two samples, no cross-sample edges
    --cat_covariates_keys batch
    --cat_covariates_embeds_nums 2
    --cat_covariates_no_edges True
    --cat_covariates_embeds_injection gene_expr_decoder

    # Training
    --lr 0.001
    --lambda_edge_recon 5000000.
    --lambda_gene_expr_recon 3000.
    --lambda_l1_masked 0.
    --lambda_l1_addon 30.
    --edge_batch_size 512
    --n_sampled_neighbors 4
    --seed 0
)
