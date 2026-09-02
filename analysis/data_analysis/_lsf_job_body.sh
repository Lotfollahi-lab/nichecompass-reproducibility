#!/bin/bash
# The body of a NicheCompass LSF job. Not submitted directly:
# submit_lsf_sanger.sh generates the #BSUB directives and runs this.
#
# Kept as a normal committed script rather than generated text, because
# generating it meant escaping backslash-newline continuations through a
# heredoc, and getting that wrong silently collapsed ´mpirun \ -np 4 \ ...´
# onto one line. A real script cannot have that class of bug.
#
# Configuration arrives through the environment, which LSF copies from the
# submitting shell. Everything after the options is forwarded to the training
# script.

set -euo pipefail

: "${N_GPUS:?N_GPUS must be exported by the submitter}"
: "${N_CORES_PER_GPU:=6}"
: "${MODEL_LABEL:?MODEL_LABEL must be exported by the submitter}"
: "${OPENMPI_MODULE:=ISG/experimental/fg12/openmpi/5.0.4-cuda12.1-lsf}"
: "${MODULES_INIT:=/usr/share/modules/init/sh}"
: "${ARGS_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# The affinity hostfile LSF writes confuses the MPI ranks, as terra also found
unset LSB_AFFINITY_HOSTFILE

echo "Host: $(hostname)"
echo "LSB_JOBID: ${LSB_JOBID:-unset}"
echo "LSB_MCPU_HOSTS: ${LSB_MCPU_HOSTS:-unset}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,memory.total,compute_mode --format=csv || true

cd "${ARGS_DIR}"

# Activate the python environment. A virtualenv (uv, venv, virtualenv) is
# preferred, since that is what this farm uses: terra and squint both keep
# theirs under /nfs/team361/sb75/.venvs/<project>. Conda is the fallback, for
# the environment.yaml in envs/.
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
    echo "Set VENV_PATH to a virtualenv, for example" >&2
    echo "  VENV_PATH=/nfs/team361/sb75/.venvs/nichecompass" >&2
    echo "or set CONDA_ENV with conda on PATH." >&2
    exit 1
fi
echo "python: $(command -v python)"

# nvidia-smi reports the host's GPUs through NVML, not this job's allocation,
# so the allocation is asserted through torch. A mismatch fails in seconds
# rather than as an opaque 'invalid device ordinal' on the other ranks after
# the data has loaded.
python - "${N_GPUS}" <<'PYEOF'
import sys
import torch
expected = int(sys.argv[1])
visible = torch.cuda.device_count()
print(f"visible CUDA devices: {visible} (expected {expected})")
if visible != expected:
    sys.exit(f"the job was allocated {visible} devices, not {expected}")
PYEOF

# The shared run configuration, also used by any other submitter
source "${ARGS_DIR}/xenium_humanppi_args.sh"

if [ "${N_GPUS}" -gt 1 ]; then
    # The module system has to be initialised explicitly, as terra does
    if [ ! -r "${MODULES_INIT}" ]; then
        echo "ERROR: cannot read ${MODULES_INIT}, so the OpenMPI module" >&2
        echo "cannot be loaded and mpirun would not be the LSF aware one." >&2
        echo "Set MODULES_INIT if it lives elsewhere on this cluster." >&2
        exit 1
    fi
    . "${MODULES_INIT}"
    module load "${OPENMPI_MODULE}"

    # NCCL settings taken from terra's parallel job. NVLS is disabled because
    # the multicast setup fails on this fabric; ring and tree collectives are
    # correct and only marginally slower.
    export NCCL_DEBUG="${NCCL_DEBUG:-WARN}"
    export NCCL_NVLS_ENABLE=0
    export NCCL_IB_DISABLE="${NCCL_IB_DISABLE:-1}"

    # mpirun does not set the rendezvous variables that torch.distributed
    # reads, so they are derived from the allocation here and forwarded to
    # every rank with -x. The port is derived from the job id so that two jobs
    # sharing a node cannot collide.
    export MASTER_ADDR="$(echo "${LSB_MCPU_HOSTS}" | awk '{print $1}')"
    export MASTER_PORT="$(( 20000 + ${LSB_JOBID:-0} % 20000 ))"
    echo "rendezvous: ${MASTER_ADDR}:${MASTER_PORT}"

    export OMP_NUM_THREADS="${N_CORES_PER_GPU}"
    export MKL_NUM_THREADS="${N_CORES_PER_GPU}"

    mpirun \
        -np "${N_GPUS}" \
        -bind-to none \
        -map-by slot \
        --mca pml ob1 --mca btl ^openib \
        -x PATH -x LD_LIBRARY_PATH -x VIRTUAL_ENV \
        -x MASTER_ADDR -x MASTER_PORT \
        -x NCCL_DEBUG -x NCCL_NVLS_ENABLE -x NCCL_IB_DISABLE \
        -x OMP_NUM_THREADS -x MKL_NUM_THREADS \
        python train_nichecompass_reference_model.py \
            "${NICHECOMPASS_ARGS[@]}" \
            --multi_gpu \
            --model_label "${MODEL_LABEL}" \
            "$@"
else
    export OMP_NUM_THREADS="$(( N_GPUS * N_CORES_PER_GPU ))"
    export MKL_NUM_THREADS="${OMP_NUM_THREADS}"
    python train_nichecompass_reference_model.py \
        "${NICHECOMPASS_ARGS[@]}" \
        --model_label "${MODEL_LABEL}" \
        "$@"
fi
