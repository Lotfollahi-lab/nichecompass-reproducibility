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

# The ´#BSUB´ lines below are read by LSF before any shell runs, so they cannot
# contain shell variables: to change the queue, the GPU count or the memory,
# either edit them here or override them on the command line, which takes
# precedence over the directives in the file:
#
#   bsub -q training-parallel -M 300000 < submit_lsf_xenium_humanppi_smoke.sh
#
# If you change the number of GPUs you have to change it in TWO places, because
# of that same restriction: the ´#BSUB -gpu "num=..."´ directive and the
# ´N_GPUS´ variable that ´torchrun´ uses. They must agree.

#BSUB -J nichecompass_xenium_smoke
#BSUB -q training-parallel
#BSUB -gpu "num=4:mode=exclusive_process:j_exclusive=yes"
#BSUB -n 16
#BSUB -R "span[hosts=1] select[mem>200000] rusage[mem=200000]"
#BSUB -M 200000
#BSUB -W 2:00
#BSUB -o logs/xenium_smoke_%J.out
#BSUB -e logs/xenium_smoke_%J.err

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

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nichecompass-reproducibility

export OMP_NUM_THREADS=$(( ${LSB_DJOB_NUMPROC:-16} / N_GPUS ))
export MKL_NUM_THREADS=${OMP_NUM_THREADS}

echo "Host: $(hostname)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv || true

# ´nvidia-smi´ reports the host's GPUs through NVML, not what this job was
# allocated, so the allocation is asserted through torch instead. A mismatch
# here fails in seconds rather than as an opaque 'invalid device ordinal' on
# ranks 1..3 after the data has loaded.
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import torch, sys; n = torch.cuda.device_count(); \
print('visible CUDA devices:', n); sys.exit(0 if n == ${N_GPUS} else 1)"

# ´--n_epochs_all_gps 0´ is the one deliberate deviation from the notebook,
# which uses 25. With only one epoch the notebook's value would leave gene
# program pruning switched off for the whole run, so the collective that keeps
# the pruning statistic identical across processes would never execute — and
# that is the single most important thing this smoke run is meant to exercise.
# Set it back to 25 for a real run.

# The same configuration the single device warm-up job uses, so the two cannot
# drift apart
source "$(dirname "$0")/xenium_humanppi_args.sh"

torchrun \
    --standalone \
    --nnodes=1 \
    --nproc_per_node=${N_GPUS} \
    train_nichecompass_reference_model.py \
    "${NICHECOMPASS_ARGS[@]}" \
    --multi_gpu \
    --model_label humanppi_1epoch_multigpu \
    --n_epochs 1 \
    --n_epochs_all_gps 0 \
    "$@"
