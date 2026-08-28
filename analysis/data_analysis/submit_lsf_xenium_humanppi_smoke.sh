#!/bin/bash
# One epoch multi-GPU smoke run of the Xenium human breast cancer configuration
# from xenium_human_breast_cancer_humanppi_gps.ipynb, on an LSF cluster.
#
#   bsub < submit_lsf_xenium_humanppi_smoke.sh
#
# The point of this job is to confirm that the distributed path works end to
# end on real data and real GPUs before committing to a long run. Every model
# and gene program argument below is pinned to the notebook, so the only
# intended differences from a full run are the epoch count and the fact that
# training is split across four devices.
#
# The data is read from ´datasets/st_data/gold/{dataset}_{batch}.h5ad´, which
# is where the notebook reads it from too.
#
# BEFORE THE FIRST RUN, populate the prior gene program caches with a SINGLE
# process run, since all four processes would otherwise race to download and
# write the same files. Running this script once without ´--multi_gpu´ and with
# ´--nproc_per_node=1´ is enough, or run the notebook's gene program cells.

#BSUB -J nichecompass_xenium_smoke
#BSUB -q gpu-normal
#BSUB -gpu "num=4:mode=exclusive_process:j_exclusive=yes"
#BSUB -n 16
#BSUB -R "span[hosts=1] select[mem>200000] rusage[mem=200000]"
#BSUB -M 200000
#BSUB -W 2:00
#BSUB -o logs/xenium_smoke_%J.out
#BSUB -e logs/xenium_smoke_%J.err

set -euo pipefail

N_GPUS=4

mkdir -p logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nichecompass-reproducibility

export OMP_NUM_THREADS=$(( ${LSB_DJOB_NUMPROC:-16} / N_GPUS ))
export MKL_NUM_THREADS=${OMP_NUM_THREADS}

echo "Host: $(hostname)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv || true

# ´--n_epochs_all_gps 0´ is the one deliberate deviation from the notebook,
# which uses 25. With only one epoch the notebook's value would leave gene
# program pruning switched off for the whole run, so the collective that keeps
# the pruning statistic identical across processes would never execute — and
# that is the single most important thing this smoke run is meant to exercise.
# Set it back to 25 for a real run.

torchrun \
    --standalone \
    --nnodes=1 \
    --nproc_per_node=${N_GPUS} \
    train_nichecompass_reference_model.py \
    --multi_gpu \
    --dataset xenium_human_breast_cancer \
    --reference_batches batch1 batch2 \
    --counts_key counts \
    --spatial_key spatial \
    --n_neighbors 8 \
    --species human \
    --model_label humanppi_1epoch_multigpu \
    `# Prior gene programs: the human PPI resource only, as in the notebook` \
    --no-include_omnipath_gps \
    --no-include_nichenet_gps \
    --no-include_mebocost_gps \
    --no-include_collectri_gps \
    --include_humanppi_gps \
    --humanppi_precision 80 \
    --humanppi_program_type intercellular \
    --humanppi_ambiguous_locality extracellular \
    --humanppi_unresolved_locality exclude \
    --humanppi_use_topology \
    --humanppi_detect_cis_complexes \
    --humanppi_min_extracellular_domain_length 30 \
    --humanppi_orient_juxtacrine_gps \
    --no-humanppi_symmetric_juxtacrine_gps \
    --humanppi_filter_ig_tcr_segments \
    --humanppi_filter_paralog_cross_pairs \
    --humanppi_min_rf_prob None \
    --humanppi_min_af_prob None \
    `# No combining, and keep only programs that still encode an interaction` \
    --gp_filter_mode none \
    --min_genes_per_gp 2 \
    --min_source_genes_per_gp 1 \
    --min_target_genes_per_gp 1 \
    `# Model architecture` \
    --conv_layer_encoder gatv2conv \
    --active_gp_thresh_ratio 0.01 \
    --n_addon_gp 100 \
    `# Batch as a categorical covariate, two samples, no cross-sample edges` \
    --cat_covariates_keys batch \
    --cat_covariates_embeds_nums 2 \
    --cat_covariates_no_edges True \
    --cat_covariates_embeds_injection gene_expr_decoder \
    `# Training: one epoch, everything else as in the notebook` \
    --n_epochs 1 \
    --n_epochs_all_gps 0 \
    --lr 0.001 \
    --lambda_edge_recon 5000000. \
    --lambda_gene_expr_recon 3000. \
    --lambda_l1_masked 0. \
    --lambda_l1_addon 30. \
    --edge_batch_size 512 \
    --n_sampled_neighbors 4 \
    --seed 0 \
    "$@"
