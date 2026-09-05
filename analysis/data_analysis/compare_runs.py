"""
Compare two NicheCompass runs of the same configuration, typically one trained on
several GPUs against one trained on a single GPU.

    python compare_runs.py \
        --run_a  .../results/humanppi_1gpu_100ep/<timestamp>/xenium_..._humanppi_1gpu_100ep.h5ad \
        --run_b  .../results/humanppi_2gpu_100ep/<timestamp>/xenium_..._humanppi_2gpu_100ep.h5ad \
        --log_a  logs/humanppi_1gpu_100ep_<jobid>.out \
        --log_b  logs/humanppi_2gpu_100ep_<jobid>.out

Two things make this comparison less obvious than it looks, and both are handled
here rather than left to the reader.

1. The latent matrix holds only the ACTIVE gene programs, in the order of
   ´adata.uns[active_gp_names_key]´ (models/nichecompass.py passes
   ´only_active_gps=True´). Two runs can prune to different sets, so the two
   matrices can have different widths and their columns do not correspond by
   position. Everything below aligns them by gene program NAME.

2. A gene program's sign is an exact gauge freedom of the training objective:
   flipping a latent dimension together with that program's decoder column
   changes nothing the model computes. Two independent runs therefore assign the
   sign independently, and roughly half the programs can come out mirrored. A
   naive correlation would read that as violent disagreement. What is compared
   here is |r|, with the number of mirrored programs reported separately -- it is
   a property of the initialisation, not a difference between the runs.
"""

import argparse
import re

import numpy as np
import pandas as pd
import anndata


def parse_log(path):
    """Pull the training time, the epoch reached and the final metrics out of a job log."""
    if path is None:
        return {}
    text = open(path, encoding="utf-8", errors="replace").read()
    out = {}
    m = re.search(r"Model training finished after (\d+) min (\d+) sec", text)
    if m:
        out["train_seconds"] = int(m.group(1)) * 60 + int(m.group(2))
    m = re.search(r"Using best model state, which was in epoch (\d+)", text)
    if m:
        out["best_epoch"] = int(m.group(1))
    epochs = re.findall(r"Epoch (\d+)/(\d+)", text)
    if epochs:
        out["epochs_run"] = int(epochs[-1][0])
        out["epochs_requested"] = int(epochs[-1][1])
    for key, label in [("auroc", "val AUROC score"),
                       ("auprc", "val AUPRC score"),
                       ("target_mse", "val target rna MSE score"),
                       ("source_mse", "val source rna MSE score")]:
        m = re.search(re.escape(label) + r":\s*([0-9.]+)", text)
        if m:
            out[key] = float(m.group(1))
    if "train_seconds" in out and out.get("epochs_run"):
        out["seconds_per_epoch"] = out["train_seconds"] / out["epochs_run"]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run_a", required=True, help="result h5ad of the first run")
    parser.add_argument("--run_b", required=True, help="result h5ad of the second run")
    parser.add_argument("--log_a", default=None, help="LSF .out log of the first run")
    parser.add_argument("--log_b", default=None, help="LSF .out log of the second run")
    parser.add_argument("--label_a", default="run A")
    parser.add_argument("--label_b", default="run B")
    parser.add_argument("--latent_key", default="nichecompass_latent")
    parser.add_argument("--active_gp_names_key", default="nichecompass_active_gp_names")
    parser.add_argument("--cluster", action="store_true",
                        help="also cluster both latents and report the ARI between them")
    parser.add_argument("--resolution", type=float, default=0.5)
    args = parser.parse_args()

    print(f"Loading {args.label_a}: {args.run_a}")
    a = anndata.read_h5ad(args.run_a)
    print(f"Loading {args.label_b}: {args.run_b}")
    b = anndata.read_h5ad(args.run_b)

    # ---------------------------------------------------------------- timing --
    log_a, log_b = parse_log(args.log_a), parse_log(args.log_b)
    if log_a or log_b:
        print("\n=== Training ===")
        rows = []
        for key, label in [("epochs_run", "epochs run"),
                           ("best_epoch", "best epoch"),
                           ("train_seconds", "training time (s)"),
                           ("seconds_per_epoch", "seconds per epoch"),
                           ("auroc", "val AUROC"),
                           ("auprc", "val AUPRC"),
                           ("target_mse", "val target RNA MSE"),
                           ("source_mse", "val source RNA MSE")]:
            if key in log_a or key in log_b:
                rows.append([label, log_a.get(key, "-"), log_b.get(key, "-")])
        print(pd.DataFrame(rows, columns=["", args.label_a, args.label_b])
              .to_string(index=False))
        if "seconds_per_epoch" in log_a and "seconds_per_epoch" in log_b:
            speedup = log_a["seconds_per_epoch"] / log_b["seconds_per_epoch"]
            print(f"\n  per-epoch speed-up of {args.label_b} over {args.label_a}: "
                  f"{speedup:.2f}x")
            print("  (compare per EPOCH, not per job: early stopping is on by "
                  "default with patience 8, so the two runs need not stop at the "
                  "same epoch)")

    # ------------------------------------------------------- active programs --
    names_a = list(a.uns[args.active_gp_names_key])
    names_b = list(b.uns[args.active_gp_names_key])
    set_a, set_b = set(names_a), set(names_b)
    shared = set_a & set_b
    jaccard = len(shared) / len(set_a | set_b) if (set_a | set_b) else float("nan")

    print("\n=== Active gene programs ===")
    print(f"  {args.label_a}: {len(set_a)}")
    print(f"  {args.label_b}: {len(set_b)}")
    print(f"  shared:  {len(shared)}   Jaccard: {jaccard:.4f}")
    only_a, only_b = sorted(set_a - set_b), sorted(set_b - set_a)
    if only_a:
        print(f"  only in {args.label_a} ({len(only_a)}): {', '.join(only_a[:8])}"
              f"{' ...' if len(only_a) > 8 else ''}")
    if only_b:
        print(f"  only in {args.label_b} ({len(only_b)}): {', '.join(only_b[:8])}"
              f"{' ...' if len(only_b) > 8 else ''}")

    # ------------------------------------------------------------- the latent --
    if not a.obs_names.equals(b.obs_names):
        common = a.obs_names.intersection(b.obs_names)
        print(f"\n  note: cell order differs; comparing the {len(common):,} shared cells")
        a, b = a[common], b[common]

    idx_a = {n: i for i, n in enumerate(names_a)}
    idx_b = {n: i for i, n in enumerate(names_b)}
    ordered = [n for n in names_a if n in shared]
    Za = np.asarray(a.obsm[args.latent_key])[:, [idx_a[n] for n in ordered]]
    Zb = np.asarray(b.obsm[args.latent_key])[:, [idx_b[n] for n in ordered]]

    # Pearson r per gene program, computed column-wise
    Za_c = Za - Za.mean(0, keepdims=True)
    Zb_c = Zb - Zb.mean(0, keepdims=True)
    denom = np.linalg.norm(Za_c, axis=0) * np.linalg.norm(Zb_c, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(denom > 0, (Za_c * Zb_c).sum(0) / denom, np.nan)
    abs_r = np.abs(r)
    mirrored = int(np.nansum(r < 0))

    print("\n=== Latent agreement, per shared gene program ===")
    print(f"  programs compared:    {len(ordered)}")
    print(f"  median |r|:           {np.nanmedian(abs_r):.4f}")
    print(f"  mean |r|:             {np.nanmean(abs_r):.4f}")
    print(f"  |r| >= 0.9:           {int(np.nansum(abs_r >= 0.9))} "
          f"({100 * np.nanmean(abs_r >= 0.9):.1f}%)")
    print(f"  |r| >= 0.7:           {int(np.nansum(abs_r >= 0.7))} "
          f"({100 * np.nanmean(abs_r >= 0.7):.1f}%)")
    print(f"  |r| <  0.3:           {int(np.nansum(abs_r < 0.3))}")
    print(f"  sign-mirrored:        {mirrored} of {len(ordered)} "
          f"({100 * mirrored / max(len(ordered), 1):.0f}%)")
    print("    ^ expected to sit near 50%. The sign is a gauge choice, not a")
    print("      disagreement: it is fixed by initialisation, and flipping it")
    print("      changes nothing the model computes.")

    worst = np.argsort(np.nan_to_num(abs_r))[:8]
    print("\n  weakest agreement:")
    for i in worst:
        print(f"    |r| = {abs_r[i]:.4f}   {ordered[i]}")

    # --------------------------------------------------------------- niches ---
    if args.cluster:
        try:
            import scanpy as sc
            from sklearn.metrics import adjusted_rand_score
            print("\n=== Niche agreement ===")
            labels = {}
            for tag, adata, Z in [(args.label_a, a, Za), (args.label_b, b, Zb)]:
                tmp = anndata.AnnData(Z)
                sc.pp.neighbors(tmp, use_rep="X")
                sc.tl.leiden(tmp, resolution=args.resolution,
                             key_added="niche", flavor="igraph", n_iterations=2)
                labels[tag] = tmp.obs["niche"].values
                print(f"  {tag}: {len(set(labels[tag]))} niches")
            ari = adjusted_rand_score(labels[args.label_a], labels[args.label_b])
            print(f"  adjusted Rand index: {ari:.4f}")
            print("    ^ clustered on the SHARED programs only, so this is the")
            print("      agreement of the niches the two runs would give you.")
        except ImportError as error:
            print(f"\n  clustering skipped: {error}")

    print("\n=== Reading this ===")
    print("  These runs are not expected to be identical. Each process draws its")
    print("  own negative edges, so a multi-GPU run differs from a single-GPU one")
    print("  the way two single-GPU runs with different seeds differ. To know")
    print("  whether the numbers above are good, run the same comparison between")
    print("  two SINGLE-GPU runs that differ only in --seed: that is the floor,")
    print("  and multi-GPU should not be meaningfully worse than it.")


if __name__ == "__main__":
    main()
