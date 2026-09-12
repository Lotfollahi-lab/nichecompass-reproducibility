"""
Populate the prior gene program caches, and nothing else.

    python prepare_gp_caches.py --gp_data_folder_path /somewhere/gp_data

Run this on a machine that HAS internet -- a laptop, a login node, the Sanger
farm -- then copy the resulting files to the cluster that does not. It needs no
spatial data, no GPU and no trained model; it only fetches the prior resources
and writes them in the exact form the training script reads back.

That last point is why this script exists rather than a note saying "download
the files". The four cached files are not the upstream downloads:

  humanppi_network_<precision>.csv   the interaction table extracted from the
                                     released archive and rewritten as TSV
                                     WITHOUT the leading '#' comment lines. The
                                     cache reader does not pass comment='#', so
                                     a raw extracted file is misparsed.
  humanppi_protein_topology.tsv      built from roughly 12,000 UniProt REST
                                     queries, with derived columns
                                     (is_gpi_anchored, the longest cytoplasmic
                                     domain, GO molecular functions). There is
                                     no upstream file that corresponds to it.
  complex_portal_human.tsv           the EBI table, saved as retrieved.
  omnipath_intercell_annotation.tsv  built from the omnipath client.

Every one is tab separated. Each is written with a .provenance.json sidecar
recording where it came from and when, which the training run reports and
validates, so copy those across too.

The arguments mirror the ones the training script takes, because the
interactome cache is keyed on the precision and the topology cache is only
complete for the accessions the chosen interaction set mentions.
"""

import argparse
import os

from nichecompass.utils import extract_gp_dict_from_humanppi_interactions


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gp_data_folder_path", type=str, required=True,
                        help="Folder to write the caches into. Created if it "
                             "does not exist.")
    parser.add_argument("--species", type=str, default="human")
    parser.add_argument("--humanppi_precision", type=str, default="80",
                        help="Must match the training run: the interactome "
                             "cache is keyed on it.")
    parser.add_argument("--humanppi_program_type", type=str,
                        default="intercellular")
    parser.add_argument("--humanppi_ambiguous_locality", type=str,
                        default="extracellular")
    parser.add_argument("--humanppi_unresolved_locality", type=str,
                        default="exclude")
    parser.add_argument("--humanppi_min_extracellular_domain_length",
                        type=int, default=30)
    args = parser.parse_args()

    folder = os.path.abspath(args.gp_data_folder_path)
    os.makedirs(folder, exist_ok=True)

    network = f"{folder}/humanppi_network_{args.humanppi_precision}.csv"
    topology = f"{folder}/humanppi_protein_topology.tsv"
    complex_portal = f"{folder}/complex_portal_human.tsv"
    omnipath_annotation = f"{folder}/omnipath_intercell_annotation.tsv"

    print(f"Writing the prior gene program caches into {folder}")
    print("This downloads a few hundred megabytes and makes roughly 120 "
          "UniProt requests, so it takes a few minutes.\n")

    # ´load_from_disk=False´ forces the fetch, ´save_to_disk=True´ writes the
    # caches. The gene program dictionary it returns is discarded: producing
    # the caches is the whole point, and the training run rebuilds the
    # dictionary from them.
    extract_gp_dict_from_humanppi_interactions(
        species=args.species,
        precision=args.humanppi_precision,
        program_type=args.humanppi_program_type,
        ambiguous_locality=args.humanppi_ambiguous_locality,
        unresolved_locality=args.humanppi_unresolved_locality,
        filter_ig_tcr_segments=True,
        filter_paralog_cross_pairs=True,
        use_topology=True,
        topology_file_path=topology,
        detect_cis_complexes=True,
        min_extracellular_domain_length=(
            args.humanppi_min_extracellular_domain_length),
        orient_juxtacrine_gps=True,
        omnipath_annotation_file_path=omnipath_annotation,
        symmetric_juxtacrine_gps=False,
        complex_portal_file_path=complex_portal,
        min_rf_prob=None,
        min_af_prob=None,
        load_from_disk=False,
        save_to_disk=True,
        ppi_network_file_path=network,
        plot_gp_gene_count_distributions=False)

    print("\nWrote:")
    missing = []
    for path in (network, topology, complex_portal, omnipath_annotation):
        if os.path.exists(path):
            size = os.path.getsize(path) / 1e6
            sidecar = " + provenance" if os.path.exists(
                f"{path}.provenance.json") else " (no provenance sidecar)"
            print(f"  {size:8.1f} MB  {os.path.basename(path)}{sidecar}")
        else:
            missing.append(path)
    for path in missing:
        print(f"  MISSING      {os.path.basename(path)}")
    if missing:
        raise SystemExit("Not every cache was written; see above.")

    print(f"\nCopy the whole folder to the cluster, then run training with")
    print(f"  GP_DATA_DIR=<destination> ...")


if __name__ == "__main__":
    main()
