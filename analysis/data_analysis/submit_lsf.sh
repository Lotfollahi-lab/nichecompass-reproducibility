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

# The ´#BSUB´ lines below are read by LSF before any shell runs, so they cannot
# contain shell variables: to change the queue, the GPU count or the memory,
# either edit them here or override them on the command line, which takes
# precedence over the directives in the file:
#
#   bsub -q training-parallel -M 300000 < submit_lsf.sh
#
# If you change the number of GPUs you have to change it in TWO places, because
# of that same restriction: the ´#BSUB -gpu "num=..."´ directive and the
# ´N_GPUS´ variable that ´torchrun´ uses. They must agree.

#BSUB -J nichecompass_train
#BSUB -q training-parallel
#BSUB -gpu "num=4:mode=exclusive_process:j_exclusive=yes"
#BSUB -n 16
#BSUB -R "span[hosts=1] select[mem>200000] rusage[mem=200000]"
#BSUB -M 200000
#BSUB -W 24:00
#BSUB -o logs/nichecompass_%J.out
#BSUB -e logs/nichecompass_%J.err

set -euo pipefail

N_GPUS=4

# LSF opens the -o / -e files before this script runs, so ´logs/´ has to exist
# BEFORE submission, not here. Create it once with ´mkdir -p logs´.

# Some parallel queues wrap every job in an MPI style launcher via JOB_STARTER,
# which starts the script once per slot. With 16 slots that would be 16 copies
# of torchrun and 64 processes contending for 4 GPUs, which fails as a spray of
# CUDA errors rather than as anything legible. Check for it and stop.
if [[ -n "${OMPI_COMM_WORLD_RANK:-}${PMI_RANK:-}${PMIX_RANK:-}${MPI_LOCALRANKID:-}" ]]; then
    echo "ERROR: this job was started by an MPI style launcher, so torchrun" >&2
    echo "would run once per slot. Check JOB_STARTER in 'bqueues -l" >&2
    echo "training-parallel' and 'bapp -l'." >&2
    exit 1
fi

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
# ´nvidia-smi´ reports the host's GPUs through NVML, not what this job was
# allocated, so the allocation is asserted through torch instead. A mismatch
# here fails in seconds rather than as an opaque 'invalid device ordinal' on
# ranks 1..3 after the data has loaded.
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import torch, sys; n = torch.cuda.device_count(); \
print('visible CUDA devices:', n); sys.exit(0 if n == ${N_GPUS} else 1)"

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
