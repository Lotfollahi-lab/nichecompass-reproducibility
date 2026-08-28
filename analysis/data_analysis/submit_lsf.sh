#!/bin/bash
# Submit a multi-GPU NicheCompass training run to an LSF cluster.
#
#   bsub < submit_lsf.sh
#
# Requests four GPUs on ONE host and starts one training process per GPU.
# Adjust the queue name, the memory and the wall clock limit to your cluster:
# queue names and whether ´-M´ is per job or per slot differ between sites.
#
# Before the first multi-GPU run, populate the prior gene program caches with a
# single process run, since all processes would otherwise race to download and
# write the same cache files. See the "Retrieve the prior gene program caches
# first" section of the multi-GPU user guide.

#BSUB -J nichecompass_train
#BSUB -q gpu-normal
#BSUB -gpu "num=4:mode=exclusive_process:j_exclusive=yes"
#BSUB -n 16
#BSUB -R "span[hosts=1] select[mem>200000] rusage[mem=200000]"
#BSUB -M 200000
#BSUB -W 24:00
#BSUB -o logs/nichecompass_%J.out
#BSUB -e logs/nichecompass_%J.err

set -euo pipefail

N_GPUS=4

mkdir -p logs

# Every process holds its own copy of the graph and of ´adata´, so the host
# memory the job needs is roughly the single device requirement times the
# number of processes. If the job dies with an out of memory error that does
# not mention CUDA, this is why: request more memory or use fewer GPUs.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nichecompass-reproducibility

# Keep the per process thread pools from oversubscribing the allocated cores
export OMP_NUM_THREADS=$(( ${LSB_DJOB_NUMPROC:-16} / N_GPUS ))
export MKL_NUM_THREADS=${OMP_NUM_THREADS}

echo "Host: $(hostname)"
echo "GPUs visible: ${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,memory.total --format=csv || true

# ´--standalone´ sets up the rendezvous on this one host and picks a free port,
# which avoids clashing with another job on a shared node. ´--multi_gpu´ tells
# NicheCompass to actually split the training; without it the script would run
# on every process independently, which is not what you want.
torchrun \
    --standalone \
    --nnodes=1 \
    --nproc_per_node=${N_GPUS} \
    train_nichecompass_reference_model.py \
    --multi_gpu \
    --dataset xenium_human_breast_cancer \
    --n_epochs 400 \
    --edge_batch_size 512 \
    "$@"
