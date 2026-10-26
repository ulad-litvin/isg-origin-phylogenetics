#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO

def read_reps(rep_file):
    """Read cluster representatives from txt file"""
    with open(rep_file) as f:
        reps = {line.strip() for line in f if line.strip()}
    return reps


def read_clusters(tsv_file):
    """Read MMseqs cluster TSV into dict: rep -> set(members)"""
    clusters = defaultdict(set)

    with open(tsv_file) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            clusters[rep].add(member)

    return clusters


def get_target_members(reps, clusters):
    """Get all members belonging to selected cluster reps"""
    selected = set()

    for rep in reps:
        if rep in clusters:
            selected.update(clusters[rep])
        else:
            print(f"Warning: rep {rep} not found in cluster file")

    return selected


def filter_fasta(fasta_file, output_file, selected_ids):
    """Write only selected sequences to output FASTA"""
    count = 0

    with open(output_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in selected_ids:
                SeqIO.write(record, out, "fasta")
                count += 1

    return count


def main():
    parser = argparse.ArgumentParser(description="Extract cluster members into FASTA")
    parser.add_argument("-r", "--reps", required=True, help="TXT file with cluster representatives")
    parser.add_argument("-t", "--tsv", required=True, help="MMseqs cluster TSV file")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file with all sequences")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")

    args = parser.parse_args()

    print("Reading representatives...")
    reps = read_reps(args.reps)

    print("Reading cluster TSV...")
    clusters = read_clusters(args.tsv)

    print("Selecting cluster members...")
    selected_ids = get_target_members(reps, clusters)

    print(f"Total selected sequences: {len(selected_ids)}")

    print("Filtering FASTA...")
    written = filter_fasta(args.fasta, args.output, selected_ids)

    print(f"Written sequences: {written}")
    print("Done.")


if __name__ == "__main__":
    main()