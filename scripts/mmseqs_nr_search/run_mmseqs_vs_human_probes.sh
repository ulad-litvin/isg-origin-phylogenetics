#!/usr/bin/env bash

set -euo pipefail

INPUT_DIR="${1:-}"

if [[ -z "$INPUT_DIR" ]]; then
    echo "Usage: $0 <input_folder>"
    exit 1
fi

for GENE_DIR in "$INPUT_DIR"/*/; do
    [[ -d "$GENE_DIR" ]] || continue

    echo "Checking $GENE_DIR"

    QUERY="${GENE_DIR}/human_probes.fasta"
    TARGET=$(find "$GENE_DIR" -maxdepth 1 -type f -name "*_ensembl_nr_hits.fasta" | head -n 1)

    if [[ ! -f "$QUERY" ]]; then
        echo "  → Skipping (human_probes.fasta not found)"
        continue
    fi

    if [[ -z "$TARGET" ]]; then
        echo "  → Skipping (*_ensembl_nr_hits.fasta not found)"
        continue
    fi

    echo "  → Running MMseqs"

    (
        cd "$GENE_DIR" || exit 1

        TARGET_BASENAME=$(basename "$TARGET")
        OUTPUT_BASENAME="${TARGET_BASENAME%.fasta}_vs_human.tsv"

        mkdir -p tmp

        mmseqs easy-search \
            human_probes.fasta \
            "$TARGET_BASENAME" \
            "$OUTPUT_BASENAME" \
            tmp \
            --threads 1 \
            --split-memory-limit 500G \
            -s 7.5 \
            --max-seqs 250000 \
            --exhaustive-search \
            --format-mode 4 \
            --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qcov,qlen,tstart,tend,tcov,tlen,evalue,bits"
    )

done