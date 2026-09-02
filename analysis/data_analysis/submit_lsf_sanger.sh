#!/bin/bash
# Submit a NicheCompass training run to the Sanger LSF farm.
#
#   bash submit_lsf_sanger.sh --n_epochs 1 --n_epochs_all_gps 0   # 4 GPU smoke run
#   bash submit_lsf_sanger.sh --n_epochs 400 --n_epochs_all_gps 25  # real run
#   N_GPUS=1 bash submit_lsf_sanger.sh --n_epochs 1               # single device
#   DRY_RUN=1 bash submit_lsf_sanger.sh --n_epochs 1              # inspect only
#
# This generates the #BSUB directives and then runs ´_lsf_job_body.sh´, which is
# a normal committed script rather than generated text. Two reasons:
#   - ´#BSUB´ lines are read before any shell runs and cannot reference
#     variables, so they have to be generated to keep the ´-gpu num=´ request
#     and the launcher's process count from drifting apart;
#   - the job body cannot be generated safely, because passing backslash
#     newline continuations through a heredoc collapses ´mpirun \ -np 4 \ ...´
#     onto one line. Directives are single lines and generate fine; the body is
#     a real script.
#
# Configuration is passed to the job through the environment, which LSF copies
# from this shell.
#
# The conventions follow the existing Sanger submitters in the ´terra´ and
# ´squint´ repositories: ´-G´ is required, the ´training-parallel´ queue is used
# with an advance reservation, and multi-GPU is launched with ´mpirun´ under the
# scheduler rather than with ´torchrun´.
#
# Verify against your own allocation:
#   bugroup -w | grep -w "$USER"       the groups you belong to
#   brsvs                              the reservations you can use
#   bqueues -l training-parallel       limits, and whether a JOB_STARTER exists

set -euo pipefail

# --- site configuration, override from the environment ----------------------
# s10396 is the group both terra and squint use on this farm. The esub rejects
# the job before queueing it if the group is wrong, and names the ones you have.
export LSF_GROUP="${LSF_GROUP:-s10396}"
export LSF_QUEUE="${LSF_QUEUE:-training-parallel}"
export LSF_RESERVATION="${LSF_RESERVATION:-lotfollahi-training-parallel}"
export N_GPUS="${N_GPUS:-4}"                      # GPUs, and ranks, per node
export N_CORES_PER_GPU="${N_CORES_PER_GPU:-6}"
export LSF_MEM_GB="${LSF_MEM_GB:-200}"            # per node
export LSF_GPU_MEM_MB="${LSF_GPU_MEM_MB:-80000}"  # gmem, as terra requests
export LSF_WALL="${LSF_WALL:-2:00}"
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

if [ "${N_GPUS}" -gt 1 ]; then
    export MODEL_LABEL="${MODEL_LABEL:-humanppi_multigpu}"
else
    export MODEL_LABEL="${MODEL_LABEL:-humanppi_singlegpu}"
fi

RESERVATION_DIRECTIVE=""
if [ -n "${LSF_RESERVATION}" ]; then
    RESERVATION_DIRECTIVE="#BSUB -U ${LSF_RESERVATION}"
fi

# Quote the forwarded arguments so that the generated command line is exactly
# what was passed, whatever it contains
FORWARDED=""
if [ "$#" -gt 0 ]; then
    FORWARDED="$(printf '%q ' "$@")"
fi

echo "Submitting to ${LSF_QUEUE} as group ${LSF_GROUP}"
echo "  GPUs/ranks  : ${N_GPUS}"
echo "  cores       : ${N_CORES} (ptile ${N_CORES})"
echo "  memory      : ${LSF_MEM_GB}G per node, gmem ${LSF_GPU_MEM_MB}MB per GPU"
echo "  reservation : ${LSF_RESERVATION:-none}"
echo "  environment : ${VENV_PATH:-${CONDA_ENV}}"
echo "  model label : ${MODEL_LABEL}"
echo "  extra args  : ${FORWARDED:-none}"

# Only the directives are generated, and every one of them is a single line
JOB=$(cat <<JOBEOF
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
exec bash ${ARGS_DIR}/_lsf_job_body.sh ${FORWARDED}
JOBEOF
)

if [ "${DRY_RUN}" != "0" ]; then
    echo "--- DRY RUN, the job that would be submitted ---"
    printf '%s\n' "${JOB}"
    echo "--- and the body it runs: ${ARGS_DIR}/_lsf_job_body.sh ---"
    exit 0
fi

printf '%s\n' "${JOB}" | bsub
