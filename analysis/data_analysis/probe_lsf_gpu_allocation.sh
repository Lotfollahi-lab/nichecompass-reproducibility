#!/bin/bash
# Probe what the Sanger LSF farm actually gives a multi-GPU job, before
# spending a real run on finding out. Takes a couple of minutes.
#
#   bash probe_lsf_gpu_allocation.sh
#   DRY_RUN=1 bash probe_lsf_gpu_allocation.sh    # print the job, submit nothing
#
# Uses the same resource directives as submit_lsf_sanger.sh, since the point is
# to learn how those particular directives are interpreted on this queue.
#
# The questions this answers, none of which can be worked out from outside the
# cluster:
#
#   1. Is the job started once, or once per slot? A queue-level JOB_STARTER
#      that wraps jobs in an MPI launcher would start it per slot.
#   2. Is ´num=4´ in the -gpu string per host or per task? With 24 slots those
#      readings differ by 4 GPUs against 96.
#   3. Does ´-M 200G´ mean what it looks like? Print the limit with its unit.
#   4. Can four MPI ranks bind to four distinct GPUs and all-reduce over NCCL?
#      That is everything the training run needs from the cluster.

set -euo pipefail

export LSF_GROUP="${LSF_GROUP:-s10396}"
export LSF_QUEUE="${LSF_QUEUE:-training-parallel}"
export LSF_RESERVATION="${LSF_RESERVATION:-lotfollahi-training-parallel}"
export N_GPUS="${N_GPUS:-4}"
export N_CORES_PER_GPU="${N_CORES_PER_GPU:-6}"
export LSF_MEM_GB="${LSF_MEM_GB:-200}"
export LSF_GPU_MEM_MB="${LSF_GPU_MEM_MB:-80000}"
# Defaults to the environment you submit from, since that is the one you
# have verified works. The repository env is only the fallback.
# A virtualenv is preferred, matching how terra and squint run on this farm.
# Point VENV_PATH elsewhere, or unset it and set CONDA_ENV, to use conda.
export VENV_PATH="${VENV_PATH:-/nfs/team361/sb75/.venvs/nichecompass}"
export CONDA_ENV="${CONDA_ENV:-${CONDA_DEFAULT_ENV:-}}"
export OPENMPI_MODULE="${OPENMPI_MODULE:-ISG/experimental/fg12/openmpi/5.0.4-cuda12.1-lsf}"
DRY_RUN="${DRY_RUN:-0}"

export ARGS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-${ARGS_DIR}/logs}"
mkdir -p "${LOG_DIR}"
N_CORES=$(( N_GPUS * N_CORES_PER_GPU ))

RESERVATION_DIRECTIVE=""
if [ -n "${LSF_RESERVATION}" ]; then
    RESERVATION_DIRECTIVE="#BSUB -U ${LSF_RESERVATION}"
fi

JOB=$(cat <<JOBEOF
#BSUB -J nichecompass_probe
#BSUB -q ${LSF_QUEUE}
#BSUB -G ${LSF_GROUP}
${RESERVATION_DIRECTIVE}
#BSUB -n ${N_CORES}
#BSUB -gpu "num=${N_GPUS}:gmem=${LSF_GPU_MEM_MB}:mode=exclusive_process:block=yes"
#BSUB -M ${LSF_MEM_GB}G
#BSUB -R "select[mem>${LSF_MEM_GB}G] rusage[mem=${LSF_MEM_GB}G] span[ptile=${N_CORES}]"
#BSUB -W 0:15
#BSUB -o ${LOG_DIR}/probe_%J.out
#BSUB -e ${LOG_DIR}/probe_%J.err
exec bash ${ARGS_DIR}/_lsf_probe_body.sh
JOBEOF
)

echo "Probing ${LSF_QUEUE} as group ${LSF_GROUP} with ${N_GPUS} GPUs"

if [ "${DRY_RUN}" != "0" ]; then
    printf '%s\n' "${JOB}"
    exit 0
fi
printf '%s\n' "${JOB}" | bsub
