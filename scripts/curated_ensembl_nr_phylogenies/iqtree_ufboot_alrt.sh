#!/bin/bash

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_fasta>"
    exit 1
fi

FASTA="$1"
THREADS=8

# filenames
base=$(basename "$FASTA" .fasta)
dir=$(dirname "$FASTA")

mafft_out="$dir/${base}.mafft.fasta"
trimal_out="$dir/${base}.trimal.fasta"

echo "Processing: $FASTA"

# -----------------------------
# Run MAFFT
# -----------------------------
echo "Running MAFFT..."
mafft --auto --reorder --thread $THREADS "$FASTA" > "$mafft_out"

# -----------------------------
# Run trimAl
# -----------------------------
echo "Running trimAl..."
trimal -in "$mafft_out" -out "$trimal_out" -gt 0.1

# -----------------------------
# Run IQ-TREE
# -----------------------------
echo "Running IQ-TREE..."
iqtree2 -nt $THREADS -s "$trimal_out" -m MFP -B 1000 -alrt 1000

echo "-----------------------------------"
echo "Pipeline completed successfully."