#!/usr/bin/env python3

# Usage: python remove_fasta_records.py <input.fasta> <remove_ids.txt> <output.fasta>

import sys
from Bio import SeqIO

def filter_fasta(fasta_file, remove_ids_file, output_file):
    # Load record IDs to remove
    with open(remove_ids_file) as f:
        remove_ids = set(line.strip() for line in f if line.strip())

    # Filter sequences
    records = (
        record for record in SeqIO.parse(fasta_file, "fasta")
        if record.id not in remove_ids
    )

    # Write filtered sequences
    SeqIO.write(records, output_file, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.fasta> <remove_ids.txt> <output.fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    remove_ids_file = sys.argv[2]
    output_file = sys.argv[3]

    filter_fasta(fasta_file, remove_ids_file, output_file)