#!/bin/bash
# Submit a NicheCompass training run to a Slurm cluster.
#
#   DATA_DIR=/path/to/h5ads bash submit_slurm.sh --n_epochs 100
#   DATA_DIR=... N_GPUS=2 bash submit_slurm.sh --n_epochs 100
#   DATA_DIR=... N_GPUS=1 bash submit_slurm.sh --n_epochs 1     # single device
#   DATA_DIR=... DRY_RUN=1 bash submit_slurm.sh --n_epochs 100  # inspect only
#
# DATA_DIR is required and is the folder holding the spatial omics files, read
# as '{DATA_DIR}/{dataset}_{batch}.h5ad'. It is passed through rather than
# baked in, so the data can live anywhere the compute nodes can see.
#
# This generates the #SBATCH directives and then runs _slurm_job_body.sh, which
# is a normal committed script rather than generated text. Two reasons, the
# same as for the LSF submitter:
#   - #SBATCH lines are read before any shell runs and cannot reference
#     variables, so they have to be generated to keep the GPU request and the
#     launcher's process count from drifting apart;
#   - the body cannot be generated safely, because passing backslash-newline
#     continuations through a heredoc collapses the launcher invocation onto
#     one line.
#
# ASKING FOR A PARTICULAR GPU MODEL. Clusters label GPU models in one of two
# ways and they are not interchangeable:
#   - a typed gres, requested as --gres=gpu:a100:N
#   - a node feature, requested as --constraint=a100
# Find out which this one uses before submitting:
#       sinfo -o '%20P %10G %40f'
# The %G column shows the gres (look for 'gpu:a100:4' rather than a bare
# 'gpu:4'); %f shows the features. Set GPU_GRES or GPU_CONSTRAINT to match. A
# type request the scheduler does not understand is not an error -- it is
# silently satisfied by whatever was free -- so the job body also asserts the
# model it actually got, through REQUIRE_GPU_MODEL, and fails if it is wrong.

set -euo pipefail

# --- site configuration, override from the environment ----------------------
export SLURM_PARTITION="${SLURM_PARTITION:-highgpu}"
export N_GPUS="${N_GPUS:-4}"                       # GPUs, and ranks, per node
export N_NODES="${N_NODES:-1}"
export N_CPUS_PER_GPU="${N_CPUS_PER_GPU:-6}"
export MEM_GB="${MEM_GB:-200}"                     # per node
export WALL="${WALL:-12:00:00}"
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"          # optional
# A100 by default, as a typed gres. If this cluster uses node features
# instead, set GPU_GRES="gpu:${N_GPUS}" and GPU_CONSTRAINT="a100".
export GPU_GRES="${GPU_GRES:-gpu:a100:${N_GPUS}}"
export GPU_CONSTRAINT="${GPU_CONSTRAINT:-}"
# Checked against nvidia-smi on the allocated node. Empty disables the check.
export REQUIRE_GPU_MODEL="${REQUIRE_GPU_MODEL:-A100}"
# Default to whichever environment is active in the submitting shell, so that
# "activate it, then submit" does the obvious thing. ´VIRTUAL_ENV´ is set by a
# virtualenv's activate script and ´CONDA_DEFAULT_ENV´ by conda's; either can
# still be overridden explicitly.
export VENV_PATH="${VENV_PATH:-${VIRTUAL_ENV:-}}"
export CONDA_ENV="${CONDA_ENV:-${CONDA_DEFAULT_ENV:-}}"
DRY_RUN="${DRY_RUN:-0}"

if [ -z "${DATA_DIR:-}" ]; then
    echo "ERROR: DATA_DIR is required. It is the folder holding the" >&2
    echo "spatial omics files, read as '{DATA_DIR}/{dataset}_{batch}.h5ad'." >&2
    echo "  DATA_DIR=/path/to/h5ads bash submit_slurm.sh --n_epochs 100" >&2
    exit 1
fi
export DATA_DIR
export ARGS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-${ARGS_DIR}/logs}"
mkdir -p "${LOG_DIR}"

N_CPUS=$(( N_GPUS * N_CPUS_PER_GPU ))

if [ "${N_GPUS}" -gt 1 ]; then
    export MODEL_LABEL="${MODEL_LABEL:-humanppi_multigpu}"
else
    export MODEL_LABEL="${MODEL_LABEL:-humanppi_singlegpu}"
fi

# Built as a list and joined with the empty entries REMOVED, because a blank
# line is not a comment and Slurm stops reading #SBATCH directives at the first
# line that is not one. An unset optional directive left as an empty string
# would therefore silently discard every directive below it.
DIRECTIVES=(
    "#SBATCH --job-name=nichecompass_${MODEL_LABEL}"
    "#SBATCH --partition=${SLURM_PARTITION}"
    "${SLURM_ACCOUNT:+#SBATCH --account=${SLURM_ACCOUNT}}"
    "#SBATCH --nodes=${N_NODES}"
    "#SBATCH --ntasks-per-node=1"
    "#SBATCH --gres=${GPU_GRES}"
    "${GPU_CONSTRAINT:+#SBATCH --constraint=${GPU_CONSTRAINT}}"
    "#SBATCH --cpus-per-task=${N_CPUS}"
    "#SBATCH --mem=${MEM_GB}G"
    "#SBATCH --time=${WALL}"
    "#SBATCH --output=${LOG_DIR}/${MODEL_LABEL}_%j.out"
    "#SBATCH --error=${LOG_DIR}/${MODEL_LABEL}_%j.err"
)

FORWARDED=""
if [ "$#" -gt 0 ]; then
    FORWARDED="$(printf '%q ' "$@")"
fi

echo "Submitting to ${SLURM_PARTITION}"
echo "  nodes       : ${N_NODES}"
echo "  GPUs/ranks  : ${N_GPUS} per node   (${GPU_GRES})"
[ -n "${GPU_CONSTRAINT}" ] && echo "  constraint  : ${GPU_CONSTRAINT}"
echo "  asserted    : ${REQUIRE_GPU_MODEL:-none}"
echo "  cores       : ${N_CPUS} per node"
echo "  memory      : ${MEM_GB}G per node"
echo "  wall clock  : ${WALL}"
echo "  data        : ${DATA_DIR}"
if [ -z "${VENV_PATH}" ] && [ -z "${CONDA_ENV}" ]; then
    echo "ERROR: no python environment. Activate the one you want the job to" >&2
    echo "use and resubmit, or set VENV_PATH (a virtualenv) or CONDA_ENV." >&2
    exit 1
fi
echo "  environment : ${VENV_PATH:-${CONDA_ENV}}"
echo "  model label : ${MODEL_LABEL}"
echo "  extra args  : ${FORWARDED:-none}"

# Only the directives are generated, and every one of them is a single line
JOB="#!/bin/bash"
for directive in "${DIRECTIVES[@]}"; do
    [ -n "${directive}" ] || continue
    JOB="${JOB}"$'\n'"${directive}"
done
JOB="${JOB}"$'\n'"exec bash ${ARGS_DIR}/_slurm_job_body.sh ${FORWARDED}"

if [ "${DRY_RUN}" != "0" ]; then
    echo "--- DRY RUN, the job that would be submitted ---"
    printf '%s\n' "${JOB}"
    echo "--- and the body it runs: ${ARGS_DIR}/_slurm_job_body.sh ---"
    exit 0
fi

printf '%s\n' "${JOB}" | sbatch --export=ALL
