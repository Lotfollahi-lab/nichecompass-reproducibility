#!/bin/bash
# The body of the LSF allocation probe. Not submitted directly:
# probe_lsf_gpu_allocation.sh generates the #BSUB directives and runs this.
#
# A real committed script rather than generated text, for the same reason as
# _lsf_job_body.sh: passing backslash newline continuations through a heredoc
# collapses the mpirun invocation onto one line.

set -euo pipefail

: "${N_GPUS:?N_GPUS must be exported by the submitter}"
: "${OPENMPI_MODULE:=ISG/experimental/fg12/openmpi/5.0.4-cuda12.1-lsf}"
: "${MODULES_INIT:=/usr/share/modules/init/sh}"

unset LSB_AFFINITY_HOSTFILE

echo "=== 1. how many times was this job body started? ==="
# One MARK line means the queue runs the job once, which is what the submitter
# assumes. Several mean a JOB_STARTER launches it once per slot, in which case
# mpirun would be started repeatedly.
echo "MARK pid=$$ ppid=$PPID host=$(hostname)"

echo
echo "=== 2. launcher environment before mpirun (should be empty) ==="
env | grep -E '^(OMPI_|PMI_|PMIX_|MPI_)' | sort || echo "  none, good"

echo
echo "=== 3. what LSF allocated ==="
env | grep -E '^(LSB_|LSF_)' | sort || true
echo "LSB_MCPU_HOSTS: ${LSB_MCPU_HOSTS:-unset}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"

echo
echo "=== 3b. does this batch job have a route to the internet? ==="
# The same downloads can succeed interactively and time out in a batch job,
# because proxy variables usually come from .bashrc, which a non-interactive
# shell does not read. This settles it rather than leaving it to be inferred.
for proxy_var in http_proxy https_proxy no_proxy HTTP_PROXY HTTPS_PROXY NO_PROXY; do
    [ -n "${!proxy_var:-}" ] && echo "  ${proxy_var}=${!proxy_var}"
done
if command -v curl >/dev/null 2>&1; then
    for host in https://omnipathdb.org/about?format=text \
                https://rest.uniprot.org/ \
                https://ftp.ebi.ac.uk/; do
        if curl -sS -m 10 -o /dev/null -w "  %{http_code} ${host}\n" "${host}"; then
            :
        else
            echo "  UNREACHABLE ${host}"
        fi
    done
else
    echo "  curl not available, skipping"
fi

echo
echo "=== 4. GPUs the host has, and what this job can see ==="
nvidia-smi --query-gpu=index,name,memory.total,compute_mode --format=csv || true

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
python - <<'PYEOF'
import torch
print("torch:", torch.__version__)
print("torch.cuda.device_count():", torch.cuda.device_count())
print("nccl available:", torch.distributed.is_nccl_available())
try:
    from nichecompass.train.distributed import detect_launcher  # noqa: F401
    print("nichecompass multi-GPU support: present")
except ImportError as error:
    print("nichecompass multi-GPU support: MISSING -", error)
PYEOF

echo
echo "=== 5. can ${N_GPUS} MPI ranks bind to distinct GPUs and all-reduce? ==="
if [ ! -r "${MODULES_INIT}" ]; then
    echo "  cannot read ${MODULES_INIT}, skipping the MPI test"
    exit 0
fi
. "${MODULES_INIT}"
module load "${OPENMPI_MODULE}"

export NCCL_DEBUG="${NCCL_DEBUG:-WARN}"
export NCCL_NVLS_ENABLE=0
export NCCL_IB_DISABLE="${NCCL_IB_DISABLE:-1}"
export MASTER_ADDR="$(echo "${LSB_MCPU_HOSTS}" | awk '{print $1}')"
export MASTER_PORT="$(( 20000 + ${LSB_JOBID:-0} % 20000 ))"
echo "rendezvous: ${MASTER_ADDR}:${MASTER_PORT}"

PROBE_PY="$(mktemp --suffix=.py)"
trap 'rm -f "${PROBE_PY}"' EXIT
cat > "${PROBE_PY}" <<'PYEOF'
# Exercises exactly the path NicheCompass uses: the launcher is detected from
# the MPI environment, the device is bound from the local rank, and a
# collective runs.
import torch

from nichecompass.train.distributed import (all_reduce_mean,
                                            cleanup_distributed,
                                            detect_launcher,
                                            get_local_rank,
                                            get_rank,
                                            get_world_size,
                                            init_distributed)

print("detect_launcher():", detect_launcher(), flush=True)
assert init_distributed(), "no process group was created"
value = torch.full((4,), float(get_rank()),
                   device=f"cuda:{get_local_rank()}")
all_reduce_mean(value)
expected = sum(range(get_world_size())) / get_world_size()
print(f"all_reduce OK rank={get_rank()} local_rank={get_local_rank()} "
      f"device={torch.cuda.current_device()} value={value[0].item()} "
      f"expected={expected}", flush=True)
cleanup_distributed()
PYEOF

mpirun \
    -np "${N_GPUS}" \
    -bind-to none \
    -map-by slot \
    --mca pml ob1 --mca btl ^openib \
    -x PATH -x LD_LIBRARY_PATH -x VIRTUAL_ENV \
    -x http_proxy -x https_proxy -x no_proxy \
    -x HTTP_PROXY -x HTTPS_PROXY -x NO_PROXY \
    -x MASTER_ADDR -x MASTER_PORT \
    -x NCCL_DEBUG -x NCCL_NVLS_ENABLE -x NCCL_IB_DISABLE \
    python "${PROBE_PY}" \
  || echo "  the ${N_GPUS} rank rendezvous FAILED, see the error above"

echo
echo "=== done. Now run: bjobs -l ${LSB_JOBID:-<jobid>} ==="
echo "  and check how num= was interpreted, and MEMLIMIT with its unit."
