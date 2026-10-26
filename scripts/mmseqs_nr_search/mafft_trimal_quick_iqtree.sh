#!/usr/bin/env bash
#
# Run MAFFT, TrimAl, and a quick IQ-TREE for each gene subfolder.
#
# Usage:
#   ./run_quick_trees.sh -t 16 -g 0.2 -m LG+F+G parent_gene_directory
#
# Options:
#   -t  Number of threads (default: 8)
#   -g  TrimAl gap threshold (-gt) (default: 0.1)
#   -m  IQ-TREE substitution model (default: LG+G)
#

set -euo pipefail


# ------------------------------------------------------------
# Default parameters
# ------------------------------------------------------------
THREADS=8
TRIMAL_GT=0.1
MODEL="LG+G"


# ------------------------------------------------------------
# Argument handling
# ------------------------------------------------------------
usage() {
    echo "Usage: $0 [-t threads] [-g trimal_gt] [-m iqtree_model] <parent_gene_directory>"
    exit 1
}

while getopts ":t:g:m:h" opt; do
    case $opt in
        t) THREADS="$OPTARG" ;;
        g) TRIMAL_GT="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        h) usage ;;
        \?)
            echo "❌ ERROR: Invalid option -$OPTARG"
            usage
            ;;
        :)
            echo "❌ ERROR: Option -$OPTARG requires an argument"
            usage
            ;;
    esac
done

shift $((OPTIND -1))

if [[ $# -ne 1 ]]; then
    echo "❌ ERROR: Missing parent gene directory"
    usage
fi

PARENT_DIR="$1"

if [[ ! -d "$PARENT_DIR" ]]; then
    echo "❌ ERROR: Directory not found: $PARENT_DIR"
    exit 1
fi


# ------------------------------------------------------------
# Start
# ------------------------------------------------------------
echo "🚀 Starting quick phylogenetic tree pipeline"
echo "📁 Parent directory : $PARENT_DIR"
echo "🧵 Threads          : $THREADS"
echo "✂️  TrimAl -gt       : $TRIMAL_GT"
echo "🧬 IQ-TREE model    : $MODEL"
echo


# ------------------------------------------------------------
# Loop through gene subfolders
# ------------------------------------------------------------
for GENE_DIR in "$PARENT_DIR"/*/; do
    [[ -d "$GENE_DIR" ]] || continue

    GENE=$(basename "$GENE_DIR")

    echo "================================================================="
    echo "🧬 Processing gene: $GENE"
    echo "================================================================="

    FASTA_FILE=$(ls "$GENE_DIR"/*_ensembl_nr_hits.fasta 2>/dev/null | head -n 1 || true)

    if [[ -z "$FASTA_FILE" ]]; then
        echo "⚠️  No *_ensembl_nr_hits.fasta found — skipping $GENE"
        echo
        continue
    fi

    echo "📥 Input FASTA: $(basename "$FASTA_FILE")"

    MAFFT_OUT="${FASTA_FILE%.fasta}.mafft.fasta"
    TRIMAL_OUT="${MAFFT_OUT%.fasta}.trimal.mafft.fasta"

    QUICK_TREE_DIR="${GENE_DIR}/quick_tree"
    PREFIX="${QUICK_TREE_DIR}/${GENE}"

    # ------------------------------------------------------------
    # MAFFT
    # ------------------------------------------------------------
    echo "🧩 Running MAFFT alignment"
    mafft --auto --reorder --thread "$THREADS" "$FASTA_FILE" > "$MAFFT_OUT"
    echo "   → MAFFT output: $(basename "$MAFFT_OUT")"

    # ------------------------------------------------------------
    # TrimAl
    # ------------------------------------------------------------
    echo "✂️  Trimming alignment with TrimAl (gap threshold: ${TRIMAL_GT})"
    trimal \
        -in "$MAFFT_OUT" \
        -out "$TRIMAL_OUT" \
        -gt "$TRIMAL_GT"
    echo "   → Trimmed MSA: $(basename "$TRIMAL_OUT")"

    # ------------------------------------------------------------
    # Prepare quick_tree directory
    # ------------------------------------------------------------
    echo "📂 Preparing quick_tree directory"
    mkdir -p "$QUICK_TREE_DIR"

    # ------------------------------------------------------------
    # IQ-TREE
    # ------------------------------------------------------------
    echo "🌳 Running IQ-TREE (fast mode, model: ${MODEL})"
    iqtree2 \
        -nt "$THREADS" \
        -s "$TRIMAL_OUT" \
        -m "$MODEL" \
        -fast \
        -pre "$PREFIX"

    echo "✅ Finished gene: $GENE"
    echo
done


# ------------------------------------------------------------
# Done
# ------------------------------------------------------------
echo "🎉 All genes processed successfully"