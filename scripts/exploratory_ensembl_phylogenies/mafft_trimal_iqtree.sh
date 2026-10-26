#!/bin/bash

# Usage: bash mafft_trimal_iqtree.sh <input_folder> <threads>

# Check inputs
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_folder> <threads>"
    exit 1
fi

INPUT_DIR=$1
THREADS=$2

# Find all fasta files in subdirectories
find "$INPUT_DIR" -type f \( -name "*.fasta" -o -name "*.fa" \) | while read -r INPUT; do
    SUBDIR=$(dirname "$INPUT")
    BASENAME=$(basename "$INPUT" .fasta)
    BASENAME=${BASENAME%.fa}  # handle .fa as well

    # Create provisional_phylogeny subfolder
    OUTDIR="$SUBDIR/provisional_phylogeny"
    mkdir -p "$OUTDIR"

    echo "Processing $INPUT ..."

    # Step 1: MAFFT alignment
    echo "Running MAFFT on $INPUT with $THREADS threads..."
    mafft --auto --reorder --thread "$THREADS" "$INPUT" > "$OUTDIR/${BASENAME}.mafft.fasta"

    # Step 2: TrimAl processing
    echo "Running trimAl..."
    trimal -in "$OUTDIR/${BASENAME}.mafft.fasta" -out "$OUTDIR/${BASENAME}.trimal.mafft.fasta" -gt 0.05

    # Step 3: IQ-TREE phylogeny
    echo "Running IQ-TREE with $THREADS threads..."
    iqtree2 -nt "$THREADS" -s "$OUTDIR/${BASENAME}.trimal.mafft.fasta" -B 1000 -pre "$OUTDIR/${BASENAME}"

    echo "Finished $INPUT!"
    echo "---------------------------------------------"
done

echo "All jobs finished successfully!"