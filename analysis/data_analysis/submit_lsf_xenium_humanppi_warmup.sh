#!/bin/bash
# Single GPU warm-up run of the Xenium human breast cancer configuration, on an
# LSF cluster.
#
#   mkdir -p logs
#   bsub < submit_lsf_xenium_humanppi_warmup.sh
#
# RUN THIS BEFORE THE MULTI-GPU JOB. It does two things that the multi-GPU job
# cannot do for itself:
#
#   1. Populates the prior gene program caches under datasets/gp_data/. Those
#      are downloaded on first use, and four processes reaching that code at
#      once would race to write the same files.
#   2. Validates the whole configuration on one device, so that a mistake in
#      the arguments or the data surfaces here rather than being mistaken for a
#      problem with the distributed path.
#
# It trains for one epoch, which is enough for both purposes.
#
# The BSUB lines are read before any shell runs and so cannot use variables.
# To change the queue or the memory without editing, override on submission:
#
#   bsub -q training-parallel -M 300000 < submit_lsf_xenium_humanppi_warmup.sh

#BSUB -J nichecompass_xenium_warmup
#BSUB -q training-parallel
#BSUB -gpu "num=1:mode=exclusive_process:j_exclusive=yes"
#BSUB -n 4
#BSUB -R "span[hosts=1] select[mem>100000] rusage[mem=100000]"
#BSUB -M 100000
#BSUB -W 2:00
#BSUB -o logs/xenium_warmup_%J.out
#BSUB -e logs/xenium_warmup_%J.err

set -euo pipefail

# LSF opens the -o / -e files before this script runs, so ´logs/´ has to exist
# BEFORE submission. Create it once with ´mkdir -p logs´.

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nichecompass-reproducibility

export OMP_NUM_THREADS=${LSB_DJOB_NUMPROC:-4}
export MKL_NUM_THREADS=${OMP_NUM_THREADS}

echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
python -c "import torch, sys; n = torch.cuda.device_count(); \
print('visible CUDA devices:', n); sys.exit(0 if n >= 1 else 1)"

# The same configuration the multi-GPU job uses, so the two cannot drift apart
source "$(dirname "$0")/xenium_humanppi_args.sh"

# No torchrun and no ´--multi_gpu´: this is deliberately the single device path,
# which is what makes it a clean baseline for the multi-GPU run to be compared
# against.
python train_nichecompass_reference_model.py \
    "${NICHECOMPASS_ARGS[@]}" \
    --model_label humanppi_1epoch_singlegpu \
    --n_epochs 1 \
    --n_epochs_all_gps 0 \
    "$@"
