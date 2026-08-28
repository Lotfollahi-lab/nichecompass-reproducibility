#!/bin/bash
# Submit a multi-GPU NicheCompass training run to a Slurm cluster.
#
#   sbatch submit_slurm.sh                 # one node, four GPUs
#   sbatch --nodes=2 submit_slurm.sh       # two nodes, four GPUs each
#
# Adjust the partition, the account, the memory and the wall clock limit to
# your cluster.
#
# Before the first multi-GPU run, populate the prior gene program caches with a
# single process run, since all processes would otherwise race to download and
# write the same cache files. See the "Retrieve the prior gene program caches
# first" section of the multi-GPU user guide.

#SBATCH --job-name=nichecompass_train
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nichecompass_%j.out
#SBATCH --error=logs/nichecompass_%j.err

set -euo pipefail

N_GPUS_PER_NODE=4

mkdir -p logs

# Every process holds its own copy of the graph and of ´adata´, so the host
# memory the job needs is roughly the single device requirement times the
# number of processes PER NODE. If the job dies with an out of memory error
# that does not mention CUDA, this is why.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nichecompass-reproducibility

# Keep the per process thread pools from oversubscribing the allocated cores.
# One task per node is requested and ´torchrun´ forks the per GPU processes
# below, so the cores of the task are shared between them.
export OMP_NUM_THREADS=$(( ${SLURM_CPUS_PER_TASK:-16} / N_GPUS_PER_NODE ))
export MKL_NUM_THREADS=${OMP_NUM_THREADS}

# Rendezvous on the first node of the allocation. The port is derived from the
# job id so that two jobs sharing a node do not collide.
MASTER_ADDR=$(scontrol show hostnames "${SLURM_JOB_NODELIST}" | head -n 1)
MASTER_PORT=$(( 20000 + SLURM_JOB_ID % 20000 ))
export MASTER_ADDR MASTER_PORT

echo "Nodes: ${SLURM_JOB_NODELIST}"
echo "Rendezvous: ${MASTER_ADDR}:${MASTER_PORT}"
nvidia-smi --query-gpu=index,name,memory.total --format=csv || true

# ´srun´ starts one task per node and each of those runs ´torchrun´, which in
# turn starts one training process per GPU on its node. The c10d rendezvous
# assigns the node ranks, so nothing has to be threaded through by hand.
# ´--multi_gpu´ tells NicheCompass to actually split the training; without it
# the script would run on every process independently, which is not what you
# want.
srun --kill-on-bad-exit=1 torchrun \
    --nnodes="${SLURM_NNODES}" \
    --nproc_per_node="${N_GPUS_PER_NODE}" \
    --rdzv_id="${SLURM_JOB_ID}" \
    --rdzv_backend=c10d \
    --rdzv_endpoint="${MASTER_ADDR}:${MASTER_PORT}" \
    train_nichecompass_reference_model.py \
    --multi_gpu \
    --dataset xenium_human_breast_cancer \
    --n_epochs 400 \
    --edge_batch_size 512 \
    "$@"
