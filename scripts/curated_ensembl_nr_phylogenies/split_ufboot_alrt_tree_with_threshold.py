#!/usr/bin/env python3

import copy
import argparse
from Bio import Phylo


def split_support(tree, mode, sh_threshold, bs_threshold):
    """
    Splits combined IQ-TREE support labels (format: "SHaLRT/UFboot")
    into two separate trees, applying dual-threshold filtering.

    For each internal node:
      - If BOTH SHaLRT >= sh_threshold AND UFboot >= bs_threshold -> keep the relevant value
      - If EITHER threshold is not met                            -> set node label to "0"

    mode = 'sh' -> write SH-aLRT values (or 0)
    mode = 'bs' -> write UFboot values  (or 0)
    """
    for clade in tree.find_clades():
        label = clade.name

        if label and "/" in label:
            try:
                sh_str, bs_str = label.split("/")
                sh = float(sh_str)
                bs = float(bs_str)

                if sh >= sh_threshold and bs >= bs_threshold:
                    clade.name = sh_str if mode == "sh" else bs_str
                else:
                    clade.name = "0"

            except ValueError:
                pass

    return tree


def build_output_name(input_path, mode):
    """
    Strips .treefile suffix (if present) and appends mode tag + .nwk.

    Example:
        gene1.trimal.fasta.treefile + mode='sh'
        -> gene1.trimal.fasta_SHaLRT.nwk
    """
    tag = "alrt" if mode == "sh" else "ufboot"
    base = input_path[: -len(".treefile")] if input_path.endswith(".treefile") else input_path
    return f"{base}.{tag}.nwk"


def main():
    parser = argparse.ArgumentParser(
        description="Split IQ-TREE SHaLRT/UFboot support values into two Newick trees."
    )
    parser.add_argument(
        "treefile",
        help="Path to the IQ-TREE output treefile (e.g. gene1.fasta.treefile)"
    )
    parser.add_argument(
        "--sh-threshold",
        type=float,
        default=80.0,
        metavar="FLOAT",
        help="Minimum SH-aLRT support value to retain (default: 80.0)"
    )
    parser.add_argument(
        "--bs-threshold",
        type=float,
        default=90.0,
        metavar="FLOAT",
        help="Minimum UFboot support value to retain (default: 90.0)"
    )
    args = parser.parse_args()

    # Validate threshold ranges
    if not (0.0 <= args.sh_threshold <= 100.0):
        parser.error("--sh-threshold must be between 0 and 100")
    if not (0.0 <= args.bs_threshold <= 100.0):
        parser.error("--bs-threshold must be between 0 and 100")

    input_tree = args.treefile
    output_sh  = build_output_name(input_tree, mode="sh")
    output_bs  = build_output_name(input_tree, mode="bs")

    tree    = Phylo.read(input_tree, "newick")
    tree_sh = split_support(copy.deepcopy(tree), mode="sh",
                            sh_threshold=args.sh_threshold,
                            bs_threshold=args.bs_threshold)
    tree_bs = split_support(copy.deepcopy(tree), mode="bs",
                            sh_threshold=args.sh_threshold,
                            bs_threshold=args.bs_threshold)

    Phylo.write(tree_sh, output_sh, "newick")
    Phylo.write(tree_bs, output_bs, "newick")

    print(f"Thresholds applied: SHaLRT >= {args.sh_threshold}  |  UFboot >= {args.bs_threshold}")
    print(f"Nodes failing either threshold are labelled 0 in both output trees.\n")
    print("Trees written:")
    print(f"  SHaLRT  -> {output_sh}")
    print(f"  UFboot  -> {output_bs}")


if __name__ == "__main__":
    main()