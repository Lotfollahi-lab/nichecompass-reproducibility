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

3. ´Add-on_<i>_GP´ is a SLOT, not a name. The add-on programs start unmasked and
   learn whatever they learn, so add-on 7 of one run has no reason to be the
   same programme as add-on 7 of another. Their names match across runs and mean
   nothing, which would silently drag the agreement statistics down. They are
   therefore reported separately from the prior programs and never pooled with
   them; to ask whether the two runs found the same add-on programmes you have to
   match them by decoder weights, not by name.
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
    Za_all = np.asarray(a.obsm[args.latent_key])
    Zb_all = np.asarray(b.obsm[args.latent_key])

    def agreement(ordered):
        """|Pearson r| per program, and how many came out sign-mirrored."""
        if not ordered:
            return np.array([]), np.array([]), 0
        Za = Za_all[:, [idx_a[n] for n in ordered]]
        Zb = Zb_all[:, [idx_b[n] for n in ordered]]
        Za_c = Za - Za.mean(0, keepdims=True)
        Zb_c = Zb - Zb.mean(0, keepdims=True)
        denom = np.linalg.norm(Za_c, axis=0) * np.linalg.norm(Zb_c, axis=0)
        # A program pruned to zero has a constant latent column, so its
        # correlation is undefined rather than bad. Left as nan and excluded
        # from every statistic below, and counted separately.
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(denom > 0, (Za_c * Zb_c).sum(0) / denom, np.nan)
        return r, np.abs(r), int(np.nansum(r < 0))

    def report(title, ordered, explain_sign):
        r, abs_r, mirrored = agreement(ordered)
        finite = np.isfinite(abs_r)
        n = int(finite.sum())
        print(f"\n=== {title} ===")
        print(f"  programs compared:    {n}"
              + (f"   ({len(ordered) - n} undefined, constant in one run)"
                 if len(ordered) != n else ""))
        if n == 0:
            return
        v = abs_r[finite]
        print(f"  median |r|:           {np.median(v):.4f}")
        print(f"  mean |r|:             {v.mean():.4f}")
        for thresh in (0.9, 0.7):
            print(f"  |r| >= {thresh}:           {int((v >= thresh).sum())} "
                  f"({100 * (v >= thresh).mean():.1f}%)")
        print(f"  |r| <  0.3:           {int((v < 0.3).sum())}")
        print(f"  sign-mirrored:        {mirrored} of {n} "
              f"({100 * mirrored / max(n, 1):.0f}%)")
        if explain_sign:
            print("    ^ expected near 50%. The sign is a gauge choice, not a")
            print("      disagreement: it is fixed by initialisation, and")
            print("      flipping it changes nothing the model computes.")
        order = np.argsort(np.where(finite, abs_r, np.inf))
        print("  weakest agreement:")
        for i in order[:6]:
            if np.isfinite(abs_r[i]):
                print(f"    |r| = {abs_r[i]:.4f}   {ordered[i]}")

    # Prior programs carry a real identity: the same name is the same
    # interaction with the same prior mask in both runs.
    prior = [n for n in names_a if n in shared and not n.startswith("Add-on")]
    addon = [n for n in names_a if n in shared and n.startswith("Add-on")]

    report("Latent agreement, prior gene programs", prior, explain_sign=True)

    if addon:
        report("Latent agreement, ADD-ON slots (not comparable by name)",
               addon, explain_sign=False)
        print("    ^ read this as a sanity check only, NOT as agreement.")
        print("      Add-on programs start unmasked and learn whatever they")
        print("      learn, so ´Add-on_7_GP´ in one run is not the same")
        print("      programme as ´Add-on_7_GP´ in the other. Low |r| here is")
        print("      expected and says nothing about the runs disagreeing. To")
        print("      compare them properly, match add-on columns between the")
        print("      two saved models by decoder weight similarity:")
        print("        _, W_a, _ = model_a.get_gp_data()")
        print("        _, W_b, _ = model_b.get_gp_data()")
        print("      normalise |W| per column and cosine-match the add-on ones.")
        ordered = prior
    else:
        ordered = prior

    # --------------------------------------------------------------- niches ---
    if args.cluster:
        try:
            import scanpy as sc
            from sklearn.metrics import adjusted_rand_score
            print("\n=== Niche agreement ===")
            labels = {}
            Za_p = Za_all[:, [idx_a[n] for n in prior]]
            Zb_p = Zb_all[:, [idx_b[n] for n in prior]]
            for tag, adata, Z in [(args.label_a, a, Za_p), (args.label_b, b, Zb_p)]:
                tmp = anndata.AnnData(Z)
                sc.pp.neighbors(tmp, use_rep="X")
                sc.tl.leiden(tmp, resolution=args.resolution,
                             key_added="niche", flavor="igraph", n_iterations=2)
                labels[tag] = tmp.obs["niche"].values
                print(f"  {tag}: {len(set(labels[tag]))} niches")
            ari = adjusted_rand_score(labels[args.label_a], labels[args.label_b])
            print(f"  adjusted Rand index: {ari:.4f}")
            print("    ^ clustered on the shared PRIOR programs only, so this")
            print("      is the agreement of the niches the two runs give you.")
            print("      Add-on slots are excluded: they are not the same")
            print("      programmes across runs, so including them would")
            print("      understate the agreement.")
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
