#!/bin/bash
# The body of a NicheCompass Slurm job. Not submitted directly:
# submit_slurm.sh generates the #SBATCH directives and runs this.
#
# Kept as a normal committed script rather than generated text, for the same
# reason as _lsf_job_body.sh: passing backslash-newline continuations through a
# heredoc collapses a multi-line launcher invocation onto one line, silently.
#
# Configuration arrives through the environment, which Slurm copies from the
# submitting shell by default. Everything after the options is forwarded to the
# training script.

set -euo pipefail

: "${N_GPUS:?N_GPUS must be exported by the submitter}"
: "${DATA_DIR:?DATA_DIR must be exported by the submitter}"
: "${MODEL_LABEL:?MODEL_LABEL must be exported by the submitter}"
: "${ARGS_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Slurm spells these differently between versions, and ´set -u´ turns a name
# that is not set into a crash rather than a fallback. Resolved once, here,
# with the alternatives tried in order, so that nothing below has to know which
# spelling this cluster uses -- and so that the script still runs outside a
# Slurm allocation, which is what makes it testable.
NODELIST="${SLURM_JOB_NODELIST:-${SLURM_NODELIST:-}}"
JOB_ID="${SLURM_JOB_ID:-${SLURM_JOBID:-0}}"
N_NODES_ALLOC="${SLURM_NNODES:-${SLURM_JOB_NUM_NODES:-1}}"
N_CPUS_ALLOC="${SLURM_CPUS_PER_TASK:-$(( N_GPUS * 6 ))}"

echo "Nodes:            ${NODELIST:-not under slurm}"
echo "SLURM job id:     ${JOB_ID}"
echo "GPUs per node:    ${N_GPUS}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,memory.total,compute_mode --format=csv || true

cd "${ARGS_DIR}"

# Activate the python environment. A virtualenv is preferred; conda is the
# fallback for the environment.yaml in envs/.
if [ -n "${VENV_PATH:-}" ] && [ -r "${VENV_PATH}/bin/activate" ]; then
    # shellcheck disable=SC1091
    source "${VENV_PATH}/bin/activate"
    echo "environment: ${VENV_PATH} (virtualenv)"
elif [ -n "${CONDA_ENV:-}" ] && command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
    echo "environment: ${CONDA_ENV} (conda)"
else
    echo "ERROR: no python environment to activate." >&2
    echo "Set VENV_PATH to a virtualenv, or CONDA_ENV with conda on PATH." >&2
    exit 1
fi
echo "python: $(command -v python)"

# The GPU model is asserted rather than assumed. Requesting a type through
# ´--gres´ or ´--constraint´ only works if the cluster actually labels its
# nodes that way, and a request the scheduler does not understand is silently
# satisfied by whatever was free. Better to fail in seconds than to discover
# afterwards that the numbers came from the wrong hardware.
if [ -n "${REQUIRE_GPU_MODEL:-}" ]; then
    if ! nvidia-smi --query-gpu=name --format=csv,noheader \
         | grep -qi -- "${REQUIRE_GPU_MODEL}"; then
        echo "ERROR: this job asked for ${REQUIRE_GPU_MODEL} but was given:" >&2
        nvidia-smi --query-gpu=index,name --format=csv,noheader >&2
        echo "The type request did not take. Check how this cluster labels" >&2
        echo "its GPUs:  sinfo -o '%20P %10G %40f'" >&2
        echo "and set GPU_GRES or GPU_CONSTRAINT in the submitter to match." >&2
        exit 1
    fi
    echo "gpu model:  ${REQUIRE_GPU_MODEL} (confirmed)"
fi

# nvidia-smi reports the host's GPUs through NVML, not this job's allocation,
# so the allocation is asserted through torch. A mismatch fails in seconds
# rather than as an opaque 'invalid device ordinal' after the data has loaded.
python - "${N_GPUS}" <<'PYEOF'
import sys
import torch
expected = int(sys.argv[1])
visible = torch.cuda.device_count()
print(f"visible CUDA devices: {visible} (expected {expected})")
if visible != expected:
    sys.exit(f"the job was allocated {visible} devices, not {expected}")
PYEOF

if [ ! -d "${DATA_DIR}" ]; then
    echo "ERROR: DATA_DIR does not exist or is not readable: ${DATA_DIR}" >&2
    exit 1
fi
echo "data:       ${DATA_DIR}"

# The prior gene program resources have to be cached BEFORE the job runs: every
# process reaches this code at the same moment and they would race to write the
# SAME cache files, which can leave a truncated file behind for later runs to
# load. Checking here turns a mid-run failure on every rank into one line up
# front.
source "${ARGS_DIR}/xenium_humanppi_args.sh"

GP_DATA_DIR="${GP_DATA_DIR:-${ARGS_DIR}/../../datasets/gp_data}"
EFFECTIVE_ARGS=("${NICHECOMPASS_ARGS[@]}" "$@")
HUMANPPI_PRECISION=""
for (( arg_index=0; arg_index<${#EFFECTIVE_ARGS[@]}; arg_index++ )); do
    if [ "${EFFECTIVE_ARGS[$arg_index]}" = "--humanppi_precision" ]; then
        HUMANPPI_PRECISION="${EFFECTIVE_ARGS[$((arg_index + 1))]}"
    fi
done
: "${HUMANPPI_PRECISION:=80}"

# Only a MULTI process run can race, so only a multi process run is refused.
# A single process run is precisely how the caches get populated, and blocking
# it would make the instruction in the error below impossible to follow.
MISSING_CACHES=""
for cache in "humanppi_network_${HUMANPPI_PRECISION}.csv" \
             "humanppi_protein_topology.tsv" \
             "complex_portal_human.tsv" "omnipath_intercell_annotation.tsv"; do
    if [ ! -r "${GP_DATA_DIR}/${cache}" ]; then
        MISSING_CACHES="${MISSING_CACHES} ${cache}"
    fi
done
if [ -n "${MISSING_CACHES}" ] && [ "${N_GPUS}" -gt 1 ]; then
    echo "ERROR: these prior gene program caches are missing from" >&2
    echo "  ${GP_DATA_DIR}" >&2
    for cache in ${MISSING_CACHES}; do echo "    ${cache}" >&2; done
    echo "Run once as a single process to populate them, then resubmit." >&2
    echo "That run is not refused, because one process cannot race itself:" >&2
    echo "  GP_DATA_DIR=${GP_DATA_DIR} \\" >&2
    echo "  DATA_DIR=${DATA_DIR} N_GPUS=1 bash submit_slurm.sh --n_epochs 1 $*" >&2
    exit 1
fi
if [ -n "${MISSING_CACHES}" ]; then
    echo "These prior gene program caches are not in ${GP_DATA_DIR} yet:"
    for cache in ${MISSING_CACHES}; do echo "    ${cache}"; done
    echo "This is a single process run, so it will fetch and write them."
    echo "That needs outbound network access from this node."
else
    echo "prior gene program caches: present in ${GP_DATA_DIR}"
fi

# torchrun sets RANK, WORLD_SIZE and LOCAL_RANK, which is one of the launchers
# NicheCompass detects, so nothing has to be threaded through by hand. The
# rendezvous is on the first node of the allocation; the port comes from the
# job id so that two jobs sharing a node cannot collide.
if [ -n "${NODELIST}" ] && command -v scontrol >/dev/null 2>&1; then
    MASTER_ADDR="$(scontrol show hostnames "${NODELIST}" | head -n 1)"
else
    # Not under Slurm, or no scontrol: a single node, so this host is it.
    MASTER_ADDR="$(hostname -s)"
fi
MASTER_PORT="$(( 20000 + JOB_ID % 20000 ))"
export MASTER_ADDR MASTER_PORT
echo "rendezvous: ${MASTER_ADDR}:${MASTER_PORT}"

# One task per node, and torchrun forks one process per GPU under it, so the
# task's cores are shared between them.
export OMP_NUM_THREADS="$(( N_CPUS_ALLOC / N_GPUS ))"
export MKL_NUM_THREADS="${OMP_NUM_THREADS}"

if [ "${N_GPUS}" -gt 1 ]; then
    srun --kill-on-bad-exit=1 torchrun \
        --nnodes="${N_NODES_ALLOC}" \
        --nproc_per_node="${N_GPUS}" \
        --rdzv_id="${JOB_ID}" \
        --rdzv_backend=c10d \
        --rdzv_endpoint="${MASTER_ADDR}:${MASTER_PORT}" \
        train_nichecompass_reference_model.py \
            "${NICHECOMPASS_ARGS[@]}" \
            --gp_data_folder_path "${GP_DATA_DIR}" \
            --multi_gpu \
            --data_folder_path "${DATA_DIR}" \
            --model_label "${MODEL_LABEL}" \
            "$@"
else
    python train_nichecompass_reference_model.py \
        "${NICHECOMPASS_ARGS[@]}" \
        --gp_data_folder_path "${GP_DATA_DIR}" \
        --data_folder_path "${DATA_DIR}" \
        --model_label "${MODEL_LABEL}" \
        "$@"
fi
