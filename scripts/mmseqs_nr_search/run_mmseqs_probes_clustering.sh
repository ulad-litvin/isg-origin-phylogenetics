#!/usr/bin/env bash

# Exit immediately if a command fails, unset variable is used, or a pipe fails
set -euo pipefail

# Check input arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_folder> <combined_output_fasta>"
    exit 1
fi

INPUT_FOLDER="$1"
COMBINED_FASTA="$2"

# Remove existing combined file if it exists
if [ -f "$COMBINED_FASTA" ]; then
    rm "$COMBINED_FASTA"
fi

# Find and process fasta files
find "$INPUT_FOLDER" -type f -name "*mmseqs_probes.fasta" | while read -r fasta_file; do
    echo "Processing: $fasta_file"

    fasta_dir=$(dirname "$fasta_file")
    base_name=$(basename "$fasta_file" .fasta)

    clustered_output="${fasta_dir}/${base_name}_clustered_50-90"
    tmp_dir="${fasta_dir}/tmp"

    mkdir -p "$tmp_dir"

    mmseqs easy-cluster "$fasta_file" "$clustered_output" "$tmp_dir" \
        --min-seq-id 0.5 \
        --cov-mode 0 -c 0.9 \
        -e 1e-4 \
        --cluster-mode 0 \
        --threads 8

    echo "Finished clustering: $fasta_file"
done

echo "Combining all *_rep_seq.fasta files into $COMBINED_FASTA ..."
find "$INPUT_FOLDER" -type f -name "*_rep_seq.fasta" -exec cat {} \; > "$COMBINED_FASTA"

echo "Done! Combined file saved as: $COMBINED_FASTA"