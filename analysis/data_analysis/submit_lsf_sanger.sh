#!/bin/bash
# Submit a NicheCompass training run to the Sanger LSF farm.
#
#   mkdir -p logs
#   bash submit_lsf_sanger.sh                       # 4 GPUs, 1 epoch smoke run
#   bash submit_lsf_sanger.sh --n_epochs 400 --n_epochs_all_gps 25   # real run
#   N_GPUS=1 bash submit_lsf_sanger.sh              # single device baseline
#   DRY_RUN=1 bash submit_lsf_sanger.sh             # print the bsub, submit nothing
#
# This is a submitter, not a job script: it builds the bsub from variables and
# pipes the job body in, so the GPU count only has to be set in one place. That
# is deliberate — ´#BSUB´ directives are read before any shell runs and cannot
# reference variables, so a self-submitting script is the only way to keep the
# ´-gpu num=´ request and the launcher's process count from drifting apart.
#
# The conventions here follow the existing Sanger submitters in the ´terra´ and
# ´squint´ repositories:
#   - ´-G´ is REQUIRED. Without it the esub rejects the job with "Sorry no
#     available user group specified for this job".
#   - the ´training-parallel´ queue is used together with an advance
#     reservation (´-U´), as terra's parallel job does.
#   - multi-GPU is launched with ´mpirun´, one rank per GPU, not with
#     ´torchrun´: that is what the farm's OpenMPI-with-LSF module supports and
#     what terra uses.
#
# Check these against your own allocation before a long run:
#   bugroup -w | grep -w "$USER"       which groups you are in
#   brsvs                              which reservations you can use
#   bqueues -l training-parallel       limits, and whether a JOB_STARTER exists

set -euo pipefail

# --- site configuration, override from the environment ----------------------
LSF_GROUP="${LSF_GROUP:-team361}"          # squint and terra also use s10396
LSF_QUEUE="${LSF_QUEUE:-training-parallel}"
LSF_RESERVATION="${LSF_RESERVATION:-lotfollahi-training-parallel}"
N_GPUS="${N_GPUS:-4}"                      # GPUs, and therefore ranks, per node
N_CORES_PER_GPU="${N_CORES_PER_GPU:-6}"
LSF_MEM_GB="${LSF_MEM_GB:-200}"            # per node
LSF_GPU_MEM_MB="${LSF_GPU_MEM_MB:-80000}"  # gmem, as terra requests
LSF_WALL="${LSF_WALL:-2:00}"
# Defaults to the environment you submit from, since that is the one you
# have verified works. The repository env is only the fallback.
CONDA_ENV="${CONDA_ENV:-${CONDA_DEFAULT_ENV:-nichecompass-reproducibility}}"
OPENMPI_MODULE="${OPENMPI_MODULE:-ISG/experimental/fg12/openmpi/5.0.4-cuda12.1-lsf}"
DRY_RUN="${DRY_RUN:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-${SCRIPT_DIR}/logs}"
mkdir -p "${LOG_DIR}"

N_CORES=$(( N_GPUS * N_CORES_PER_GPU ))

# Only the multi-GPU case needs MPI. One GPU runs the single device path, which
# is the same code NicheCompass has always run.
if [ "${N_GPUS}" -gt 1 ]; then
    MULTI_GPU_FLAG="--multi_gpu"
    MODEL_LABEL="${MODEL_LABEL:-humanppi_multigpu}"
else
    MULTI_GPU_FLAG=""
    MODEL_LABEL="${MODEL_LABEL:-humanppi_singlegpu}"
fi

RESERVATION_DIRECTIVE=""
if [ -n "${LSF_RESERVATION}" ]; then
    RESERVATION_DIRECTIVE="#BSUB -U ${LSF_RESERVATION}"
fi

echo "Submitting to ${LSF_QUEUE} as group ${LSF_GROUP}"
echo "  GPUs/ranks : ${N_GPUS}"
echo "  cores      : ${N_CORES} (ptile ${N_CORES})"
echo "  memory     : ${LSF_MEM_GB}G per node, gmem ${LSF_GPU_MEM_MB}MB per GPU"
echo "  reservation: ${LSF_RESERVATION:-none}"
echo "  extra args : $*"

JOB_BODY=$(cat <<JOBEOF
#BSUB -J nichecompass_${MODEL_LABEL}
#BSUB -q ${LSF_QUEUE}
#BSUB -G ${LSF_GROUP}
${RESERVATION_DIRECTIVE}
#BSUB -n ${N_CORES}
#BSUB -gpu "num=${N_GPUS}:gmem=${LSF_GPU_MEM_MB}:mode=exclusive_process:block=yes"
#BSUB -M ${LSF_MEM_GB}G
#BSUB -R "select[mem>${LSF_MEM_GB}G] rusage[mem=${LSF_MEM_GB}G] span[ptile=${N_CORES}]"
#BSUB -W ${LSF_WALL}
#BSUB -o ${LOG_DIR}/${MODEL_LABEL}_%J.out
#BSUB -e ${LOG_DIR}/${MODEL_LABEL}_%J.err

set -euo pipefail

# terra unsets this: the affinity hostfile LSF writes confuses the MPI ranks
unset LSB_AFFINITY_HOSTFILE

echo "Host: \$(hostname)"
echo "LSB_JOBID: \${LSB_JOBID}"
echo "LSB_MCPU_HOSTS: \${LSB_MCPU_HOSTS:-unset}"
echo "CUDA_VISIBLE_DEVICES: \${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,memory.total,compute_mode --format=csv || true

cd ${SCRIPT_DIR}

source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${CONDA_ENV}

# nvidia-smi reports the host's GPUs through NVML, not this job's allocation,
# so the allocation is asserted through torch. A mismatch fails in seconds
# rather than as an opaque 'invalid device ordinal' on the other ranks after
# the data has loaded.
python -c "import torch, sys; n = torch.cuda.device_count(); \\
print('visible CUDA devices:', n); sys.exit(0 if n == ${N_GPUS} else 1)"

source ${SCRIPT_DIR}/xenium_humanppi_args.sh

if [ "${N_GPUS}" -gt 1 ]; then
    . /usr/share/modules/init/sh
    module load ${OPENMPI_MODULE}

    # NCCL settings taken from terra's parallel job. NVLS is disabled because
    # the multicast setup fails on this fabric; ring and tree collectives are
    # correct and only marginally slower.
    export NCCL_DEBUG=WARN
    export NCCL_NVLS_ENABLE=0
    export NCCL_IB_DISABLE=\${NCCL_IB_DISABLE:-1}

    # mpirun does not set the rendezvous variables that torch.distributed
    # reads, so they are derived from the allocation here and forwarded to
    # every rank with -x. The port is derived from the job id so that two jobs
    # sharing a node cannot collide.
    export MASTER_ADDR=\$(echo \${LSB_MCPU_HOSTS} | awk '{print \$1}')
    export MASTER_PORT=\$(( 20000 + LSB_JOBID % 20000 ))
    echo "rendezvous: \${MASTER_ADDR}:\${MASTER_PORT}"

    export OMP_NUM_THREADS=${N_CORES_PER_GPU}
    export MKL_NUM_THREADS=${N_CORES_PER_GPU}

    mpirun \\
        -np ${N_GPUS} \\
        -bind-to none \\
        -map-by slot \\
        --mca pml ob1 --mca btl ^openib \\
        -x PATH -x LD_LIBRARY_PATH \\
        -x MASTER_ADDR -x MASTER_PORT \\
        -x NCCL_DEBUG -x NCCL_NVLS_ENABLE -x NCCL_IB_DISABLE \\
        -x OMP_NUM_THREADS -x MKL_NUM_THREADS \\
        python train_nichecompass_reference_model.py \\
            "\${NICHECOMPASS_ARGS[@]}" \\
            ${MULTI_GPU_FLAG} \\
            --model_label ${MODEL_LABEL} \\
            $*
else
    export OMP_NUM_THREADS=${N_CORES}
    export MKL_NUM_THREADS=${N_CORES}
    python train_nichecompass_reference_model.py \\
        "\${NICHECOMPASS_ARGS[@]}" \\
        --model_label ${MODEL_LABEL} \\
        $*
fi
JOBEOF
)

if [ "${DRY_RUN}" != "0" ]; then
    echo "--- DRY RUN, the job that would be submitted ---"
    printf '%s\n' "${JOB_BODY}"
    exit 0
fi

printf '%s\n' "${JOB_BODY}" | bsub
