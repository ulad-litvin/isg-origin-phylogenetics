#!/bin/bash

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <curated_aa.fasta> <curated_nt.fasta>"
    exit 1
fi

AA_FASTA="$1"
NT_FASTA="$2"
THREADS=8

# Filenames
aa_base=$(basename "$AA_FASTA" .fasta)
nt_base=$(basename "$NT_FASTA" .fasta)
dir=$(dirname "$AA_FASTA")

# Output files
mafft_aa="$dir/${aa_base}.mafft.fasta"
pal2nal_out="$dir/${aa_base}.pal2nal.fasta"
trimal_out="$dir/${aa_base}.pal2nal.trimal.fasta"

echo "Processing AA FASTA: $AA_FASTA"
echo "Processing NT FASTA: $NT_FASTA"

# -----------------------------
# Run MAFFT on amino acid sequences
# -----------------------------
echo "Running MAFFT on amino acid sequences..."
mafft --auto --reorder --thread $THREADS "$AA_FASTA" > "$mafft_aa"

# -----------------------------
# Run PAL2NAL
# -----------------------------
echo "Running PAL2NAL to convert AA alignment to codon alignment..."
pal2nal.pl "$mafft_aa" "$NT_FASTA" -output fasta > "$pal2nal_out"

# -----------------------------
# Run trimAl on codon alignment
# -----------------------------
echo "Running trimAl on codon alignment..."
trimal -in "$pal2nal_out" -out "$trimal_out" -gt 0.1

# -----------------------------
# Run IQ-TREE on trimmed codon alignment
# -----------------------------
echo "Running IQ-TREE (fast nucleotide tree)..."
iqtree2 -nt $THREADS -s "$trimal_out" -m GTR+G -fast

echo "-----------------------------------"
echo "Pipeline completed successfully."