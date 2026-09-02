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

LSF_GROUP="${LSF_GROUP:-team361}"
LSF_QUEUE="${LSF_QUEUE:-training-parallel}"
LSF_RESERVATION="${LSF_RESERVATION:-lotfollahi-training-parallel}"
N_GPUS="${N_GPUS:-4}"
N_CORES_PER_GPU="${N_CORES_PER_GPU:-6}"
LSF_MEM_GB="${LSF_MEM_GB:-200}"
LSF_GPU_MEM_MB="${LSF_GPU_MEM_MB:-80000}"
# Defaults to the environment you submit from, since that is the one you
# have verified works. The repository env is only the fallback.
CONDA_ENV="${CONDA_ENV:-${CONDA_DEFAULT_ENV:-nichecompass-reproducibility}}"
OPENMPI_MODULE="${OPENMPI_MODULE:-ISG/experimental/fg12/openmpi/5.0.4-cuda12.1-lsf}"
DRY_RUN="${DRY_RUN:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-${SCRIPT_DIR}/logs}"
mkdir -p "${LOG_DIR}"
N_CORES=$(( N_GPUS * N_CORES_PER_GPU ))

RESERVATION_DIRECTIVE=""
if [ -n "${LSF_RESERVATION}" ]; then
    RESERVATION_DIRECTIVE="#BSUB -U ${LSF_RESERVATION}"
fi

JOB_BODY=$(cat <<JOBEOF
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

unset LSB_AFFINITY_HOSTFILE

echo "=== 1. how many times was this job body started? ==="
# One MARK line means the queue runs the job once, which is what the submitter
# assumes. Several mean a JOB_STARTER launches it per slot.
echo "MARK pid=\$\$ ppid=\$PPID host=\$(hostname)"

echo
echo "=== 2. launcher environment before mpirun (should be empty) ==="
env | grep -E '^(OMPI_|PMI_|PMIX_|MPI_)' | sort || echo "  none, good"

echo
echo "=== 3. what LSF allocated ==="
env | grep -E '^(LSB_|LSF_)' | sort
echo "LSB_MCPU_HOSTS: \${LSB_MCPU_HOSTS:-unset}"
echo "CUDA_VISIBLE_DEVICES: \${CUDA_VISIBLE_DEVICES:-unset}"

echo
echo "=== 4. GPUs the host has, and what this job can see ==="
nvidia-smi --query-gpu=index,name,memory.total,compute_mode --format=csv || true

source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${CONDA_ENV}
python -c "
import torch
print('torch:', torch.__version__)
print('torch.cuda.device_count():', torch.cuda.device_count())
print('nccl available:', torch.distributed.is_nccl_available())
"

echo
echo "=== 5. can ${N_GPUS} MPI ranks bind to distinct GPUs and all-reduce? ==="
. /usr/share/modules/init/sh
module load ${OPENMPI_MODULE}

export NCCL_DEBUG=WARN
export NCCL_NVLS_ENABLE=0
export NCCL_IB_DISABLE=\${NCCL_IB_DISABLE:-1}
export MASTER_ADDR=\$(echo \${LSB_MCPU_HOSTS} | awk '{print \$1}')
export MASTER_PORT=\$(( 20000 + LSB_JOBID % 20000 ))
echo "rendezvous: \${MASTER_ADDR}:\${MASTER_PORT}"

cat > /tmp/probe_ddp_\${LSB_JOBID}.py <<'PYEOF'
# Exercises exactly the path NicheCompass uses: the launcher is detected from
# the MPI environment, the device is bound from the local rank, and a
# collective is run.
from nichecompass.train.distributed import (detect_launcher,
                                            get_local_rank,
                                            get_rank,
                                            get_world_size,
                                            init_distributed,
                                            all_reduce_mean,
                                            cleanup_distributed)
import torch

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

mpirun \\
    -np ${N_GPUS} \\
    -bind-to none -map-by slot \\
    --mca pml ob1 --mca btl ^openib \\
    -x PATH -x LD_LIBRARY_PATH \\
    -x MASTER_ADDR -x MASTER_PORT \\
    -x NCCL_DEBUG -x NCCL_NVLS_ENABLE -x NCCL_IB_DISABLE \\
    python /tmp/probe_ddp_\${LSB_JOBID}.py \\
  || echo "  the ${N_GPUS} rank rendezvous FAILED, see the error above"
rm -f /tmp/probe_ddp_\${LSB_JOBID}.py

echo
echo "=== done. Now run: bjobs -l \${LSB_JOBID} ==="
echo "  and check how num= was interpreted and what MEMLIMIT is, with its unit."
JOBEOF
)

if [ "${DRY_RUN}" != "0" ]; then
    printf '%s\n' "${JOB_BODY}"
    exit 0
fi
printf '%s\n' "${JOB_BODY}" | bsub
