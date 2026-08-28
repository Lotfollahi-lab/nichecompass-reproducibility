#!/bin/bash
# Probe what an LSF queue actually gives a multi-GPU job, before spending a real
# run on finding out. Takes seconds of GPU time.
#
#   mkdir -p logs
#   bsub < probe_lsf_gpu_allocation.sh
#   # then, once it has run:
#   bjobs -l -gpu <jobid>
#
# Keep the resource directives here IDENTICAL to the real job, since the point
# is to learn how these particular directives are interpreted on this queue.
#
# The four questions this answers, none of which can be worked out from
# outside the cluster:
#
#   1. Does the queue wrap jobs in an MPI style launcher? Exactly one MARK line
#      means no. Sixteen MARK lines mean the script is started once per slot,
#      and torchrun must NOT be used as written, because it would start 4
#      processes per copy.
#   2. Is ´num=4´ in the -gpu string per host or per task? With -n 16 those
#      readings differ by 4 GPUs against 64.
#   3. Is the memory limit in KB or MB? -M 200000 is either 200 GB or 195 MB,
#      and the second kills the job seconds in, looking like an application bug.
#   4. Does the job land on one host, and does torch see the GPUs?

#BSUB -J nichecompass_probe
#BSUB -q training-parallel
#BSUB -gpu "num=4:mode=exclusive_process:j_exclusive=yes"
#BSUB -n 16
#BSUB -R "span[hosts=1] select[mem>200000] rusage[mem=200000]"
#BSUB -M 200000
#BSUB -W 0:10
#BSUB -o logs/probe_%J.out
#BSUB -e logs/probe_%J.err

echo "=== 1. how many times was this script started? ==="
# One line means the queue runs the job once. Several mean a JOB_STARTER is
# launching it per task, which torchrun cannot coexist with.
echo "MARK pid=$$ ppid=$PPID host=$(hostname)"
ps -o pid,ppid,args -p $PPID 2>/dev/null || true

echo
echo "=== 2. launcher environment (any of these means a launcher is in play) ==="
env | grep -E '^(OMPI_|PMI_|PMIX_|MPI_)' | sort || echo "  none, good"

echo
echo "=== 3. what LSF allocated ==="
env | grep -E '^(LSB_|LSF_)' | sort

echo
echo "=== 4. GPUs: what the host has, and what this job can see ==="
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,memory.total,compute_mode --format=csv || true

source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate nichecompass-reproducibility 2>/dev/null || true
# torch.cuda.device_count() reflects the allocation; nvidia-smi reflects the
# host, so these two disagreeing is itself the answer to question 2
python -c "
import torch
print('torch:', torch.__version__)
print('torch.cuda.is_available():', torch.cuda.is_available())
print('torch.cuda.device_count():', torch.cuda.device_count())
print('nccl available:', torch.distributed.is_nccl_available())
" || true

echo
echo "=== 5. can four processes actually rendezvous and reduce? ==="
# The smallest end to end test of the collective path: if this prints
# 'all_reduce OK' four times, NCCL, the device binding and the rendezvous all
# work, which is everything the training run needs from the cluster.
cat > /tmp/probe_ddp_$$.py <<'PYEOF'
import os
import torch
import torch.distributed as dist

local_rank = int(os.environ["LOCAL_RANK"])
torch.cuda.set_device(local_rank)
dist.init_process_group("nccl")
tensor = torch.full((4,), float(dist.get_rank()), device=f"cuda:{local_rank}")
dist.all_reduce(tensor, op=dist.ReduceOp.SUM)
expected = sum(range(dist.get_world_size()))
print(f"all_reduce OK rank={dist.get_rank()} local_rank={local_rank} "
      f"device={torch.cuda.current_device()} value={tensor[0].item()} "
      f"expected={expected}")
dist.barrier(device_ids=[local_rank])
dist.destroy_process_group()
PYEOF
torchrun --standalone --nnodes=1 --nproc_per_node=4 /tmp/probe_ddp_$$.py || \
    echo "  the four process rendezvous FAILED, see the error above"
rm -f /tmp/probe_ddp_$$.py

echo
echo "=== done. Now run: bjobs -l -gpu \$LSB_JOBID ==="
echo "  and check how num= was interpreted and what MEMLIMIT is, with its unit."
