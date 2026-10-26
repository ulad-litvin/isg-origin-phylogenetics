#!/bin/bash

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_fasta> <min_seq_id> <coverage>"
    exit 1
fi

FASTA="$1"
MIN_ID="$2"
COV="$3"

THREADS=8
TMP_DIR="tmp"

gene_dir=$(dirname "$FASTA")
filename=$(basename "$FASTA" .fasta)

outdir="$gene_dir/mmseqs_clusters"
fasttree_dir="$gene_dir/fast_tree"

mkdir -p "$outdir"
mkdir -p "$fasttree_dir"

base_name="${filename%_processed}"
level="${MIN_ID}-${COV}"

echo "Processing: $FASTA"
echo "MMseqs parameters: min_seq_id=$MIN_ID, coverage=$COV"

# -----------------------------
# Run MMseqs clustering
# -----------------------------

out="$outdir/${base_name}_${level}"

echo "Running MMseqs..."
mmseqs easy-cluster "$FASTA" "$out" "$TMP_DIR" \
    --min-seq-id "$MIN_ID" --cov-mode 1 -c "$COV" --threads $THREADS

rep_file="${out}_rep_seq.fasta"

if [ ! -f "$rep_file" ]; then
    echo "Error: clustering failed, no rep_seq file found"
    exit 1
fi

# Count representatives
rep_count=$(grep -c "^>" "$rep_file")
echo "Cluster representatives: $rep_count"

# -----------------------------
# Prepare filenames
# -----------------------------

suffix=$(basename "$rep_file" | sed -E 's/^[0-9]+_//; s/_rep_seq\.fasta$//')
base="${rep_count}_${suffix}"

mafft_out="$fasttree_dir/${base}.mafft.fasta"
trimal_out="$fasttree_dir/${base}.trimal.mafft.fasta"

# -----------------------------
# Run MAFFT
# -----------------------------

echo "Running MAFFT..."
mafft --auto --reorder --thread $THREADS "$rep_file" > "$mafft_out"

# -----------------------------
# Run trimAl
# -----------------------------

echo "Running trimAl..."
trimal -in "$mafft_out" -out "$trimal_out" -gt 0.1

# -----------------------------
# Run IQ-TREE
# -----------------------------

echo "Running IQ-TREE..."
iqtree2 -nt $THREADS -s "$trimal_out" -m LG+G -fast

echo "-----------------------------------"
echo "Pipeline completed successfully."