"""
Make an .h5ad written by anndata >= 0.12 readable by anndata < 0.12.

anndata 0.12 serialises None as a dataset with encoding-type 'null'. Older
anndata has no reader registered for it and dies on read. This deletes every
such dataset. The usual culprit is uns['log1p']['base'], which scanpy sets to
None and which nothing in NicheCompass reads.

Usage:
    python strip_nulls.py FILE.h5ad          # list what would be removed
    python strip_nulls.py FILE.h5ad --apply  # remove it, in place
"""
import shutil
import sys

import h5py


def find_nulls(path):
    found = []
    with h5py.File(path, "r") as f:
        # Collect first, mutate later: visititems must not run while the tree
        # is being changed.
        f.visititems(
            lambda name, obj: found.append(name)
            if obj.attrs.get("encoding-type") == "null" else None)
    return found


def strip(path, apply):
    nulls = find_nulls(path)
    if not nulls:
        print(f"No null-encoded values in {path} - nothing to do.")
        return 0
    print(f"Null-encoded values in {path}:")
    for name in nulls:
        print(f"  /{name}")
    if not apply:
        print("\nDry run. Re-run with --apply to remove them.")
        return 0

    backup = path + ".bak"
    print(f"\nBacking up to {backup}")
    shutil.copy2(path, backup)

    with h5py.File(path, "a") as f:
        for name in nulls:
            del f[name]
        # scanpy < 1.10 reads uns['log1p']['base'] without a .get() guard, so
        # an empty log1p dict raises KeyError where a missing one does not.
        if "uns/log1p" in f and len(f["uns/log1p"]) == 0:
            print("  removing now-empty /uns/log1p")
            del f["uns/log1p"]
    print(f"Removed {len(nulls)} null value(s).")
    print("Note: HDF5 does not reclaim freed space, so the file size will not "
          "shrink.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(strip(sys.argv[1], "--apply" in sys.argv))
