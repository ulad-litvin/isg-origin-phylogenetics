#!/bin/bash

set -e

# -----------------------------
# Usage check
# -----------------------------
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_folder>"
    echo "  input_folder: directory containing subfolders with aligned/trimmed .fasta files"
    exit 1
fi

INPUT_DIR="$1"
THREADS=8

# Validate that the input is a directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: '$INPUT_DIR' is not a directory."
    exit 1
fi

echo "Scanning for FASTA files in subfolders of: $INPUT_DIR"
echo "-----------------------------------"

# -----------------------------
# Find and process FASTA files
# -----------------------------
found=0

# Loop over every .fasta file exactly one level deep (subfolder/*.fasta)
for FASTA in "$INPUT_DIR"/*/*.fasta; do

    # Guard against the glob returning no matches
    [ -e "$FASTA" ] || { echo "No .fasta files found in any subfolder."; exit 1; }

    found=$((found + 1))
    base=$(basename "$FASTA" .fasta)
    dir=$(dirname "$FASTA")

    echo "[$found] Processing: $FASTA"

    # -----------------------------
    # Run IQ-TREE
    # -----------------------------
    echo "  Running IQ-TREE..."
    iqtree2 \
        -nt "$THREADS" \
        -s  "$FASTA" \
        -m  MFP \
        -B  1000 \
        -alrt 1000 \
        --prefix "$dir/$base"

    echo "  Done → output in: $dir/"
    echo "-----------------------------------"

done

if [ "$found" -eq 0 ]; then
    echo "Warning: no .fasta files were found under '$INPUT_DIR'."
    exit 1
fi

echo "All $found file(s) processed successfully."