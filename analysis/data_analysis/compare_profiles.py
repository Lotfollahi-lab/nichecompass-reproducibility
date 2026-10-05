"""
Compare the training speed of runs made with --profile step, model part first.

    python compare_profiles.py 1gpu=logs/speed_1gpu_<jobid>.out \
        2gpu=logs/speed_2gpu_<jobid>.out 4gpu=logs/speed_4gpu_<jobid>.out

Each argument is label=log file; the first run is the baseline. The script reads
the [nc-prof-json] line the trainer prints at the end of a profiled run and
splits every epoch into its parts, timed on the slowest process:

    model       forward passes and loss, backward with the gradient all-reduce,
                and the optimizer step (train.forward, train.backward,
                train.optimizer)
    to device   host-to-device copies of the batches (train.to_device)
    sampling    neighbour sampling and batch assembly on the CPU, training and
                validation (loader.sample)
    validation  validation forward passes, gathering and metrics
    other       the rest of the training loop: loss .item() calls, epoch
                reductions, early stopping, and anything not instrumented

The speed-up of each part is the baseline's time per epoch divided by the run's.
Compare only runs of the same configuration on the same GPU model, with the same
number of epochs, early stopping and the best-model reload switched off
(--no-use_early_stopping --no-reload_best_model) and, for a pure speed
comparison, no pruning in the timed epochs (--n_epochs_all_gps above
--n_epochs). --profile step synchronises the GPU around every timed stage, so its
epochs are slower than unprofiled ones; compare step runs only with step runs.
"""
import json
import sys

PARTS = {
    "model": ("train.forward", "train.backward", "train.optimizer"),
    "to device": ("train.to_device",),
    "sampling": ("loader.sample",),
    "validation": ("val.forward", "val.accumulate", "val.gather", "val.metrics"),
}


def read_profile(path):
    payload = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("[nc-prof-json]"):
                candidate = json.loads(line[len("[nc-prof-json]"):])
                if "probes" in candidate:
                    payload = candidate
    if payload is None:
        sys.exit(f"{path}: no [nc-prof-json] line with probes; was the run "
                 "made with --profile step?")
    if payload["mode"] != "step":
        sys.exit(f"{path}: profiled at mode={payload['mode']!r}; the per-step "
                 "stages are timed only with --profile step")
    return payload


def per_epoch(payload):
    epochs = max(payload["epochs"], 1)
    probes = payload["probes"]
    parts = {name: sum(probes.get(p, {}).get("max", 0.0) for p in members)
             / epochs for name, members in PARTS.items()}
    total = (payload["training_time_s"] or 0.0) / epochs
    parts["other"] = max(total - sum(parts.values()), 0.0)
    parts["training loop"] = total
    steps = probes.get("train.forward", {}).get("calls", 0) / epochs
    return parts, steps


def main():
    runs = []
    for arg in sys.argv[1:]:
        label, _, path = arg.partition("=")
        if not path:
            sys.exit("arguments are label=log_file")
        payload = read_profile(path)
        parts, steps = per_epoch(payload)
        runs.append((label, payload, parts, steps))
    if not runs:
        sys.exit(__doc__)

    names = list(PARTS) + ["other", "training loop"]
    print(f"{'seconds per epoch':<18}" + "".join(f"{r[0]:>14}" for r in runs))
    for key, value in (("processes", lambda r: r[1]["world_size"]),
                       ("batch convention",
                        lambda r: r[1]["meta"].get("batch_size_scaling", "?")),
                       ("epochs", lambda r: r[1]["epochs"]),
                       ("steps per epoch", lambda r: f"{r[3]:.0f}")):
        print(f"{key:<18}" + "".join(f"{str(value(r)):>14}" for r in runs))
    for name in names:
        print(f"{name:<18}" + "".join(f"{r[2][name]:>14.2f}" for r in runs))
    base = runs[0][2]
    print(f"\nspeed-up against {runs[0][0]}")
    for name in names:
        cells = []
        for r in runs:
            cells.append(f"{base[name] / r[2][name]:>13.2f}x" if r[2][name] > 0
                         else f"{'–':>14}")
        print(f"{name:<18}" + "".join(cells))
    print("\nmodel time per optimizer step (ms)")
    print(f"{'':<18}" + "".join(
        f"{1000 * r[2]['model'] / r[3]:>14.1f}" if r[3] else f"{'–':>14}"
        for r in runs))


if __name__ == "__main__":
    main()
