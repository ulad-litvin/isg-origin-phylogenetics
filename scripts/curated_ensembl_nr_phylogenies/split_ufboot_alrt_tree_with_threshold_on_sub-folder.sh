#!/bin/bash

set -e

# -----------------------------
# Usage check
# -----------------------------
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_folder> [--sh-threshold FLOAT] [--bs-threshold FLOAT]"
    echo ""
    echo "  input_folder          Directory to search recursively for .treefile files"
    echo "  --sh-threshold FLOAT  Minimum SH-aLRT support value (default: 80)"
    echo "  --bs-threshold FLOAT  Minimum UFboot support value  (default: 90)"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/trees --sh-threshold 85 --bs-threshold 95"
    exit 1
fi

PYTHON_SCRIPT="/home2/2820395l/documents/dai_evolution_phylogenies/scripts/final_trees/split_ufboot_alrt_tree_with_threshold.py"

# -----------------------------
# Defaults
# -----------------------------
SH_THRESHOLD=80
BS_THRESHOLD=90

# -----------------------------
# Parse arguments
# -----------------------------
INPUT_DIR="$1"
shift  # move past the positional argument so only flags remain in "$@"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sh-threshold)
            SH_THRESHOLD="$2"
            shift 2
            ;;
        --bs-threshold)
            BS_THRESHOLD="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown argument '$1'"
            echo "Usage: $0 <input_folder> [--sh-threshold FLOAT] [--bs-threshold FLOAT]"
            exit 1
            ;;
    esac
done

# -----------------------------
# Validate inputs
# -----------------------------
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: '$INPUT_DIR' is not a directory."
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: '$PYTHON_SCRIPT' not found in the current directory."
    echo "Make sure split_tree.py is in the same folder as this script."
    exit 1
fi

# Ensure thresholds are numeric
if ! [[ "$SH_THRESHOLD" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --sh-threshold must be a positive number (got: '$SH_THRESHOLD')"
    exit 1
fi

if ! [[ "$BS_THRESHOLD" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: --bs-threshold must be a positive number (got: '$BS_THRESHOLD')"
    exit 1
fi

echo "Searching for .treefile files in: $INPUT_DIR"
echo "SH-aLRT threshold : $SH_THRESHOLD"
echo "UFboot threshold   : $BS_THRESHOLD"
echo "-----------------------------------"

# -----------------------------
# Find and process treefiles
# -----------------------------
found=0

while IFS= read -r -d '' TREEFILE; do

    found=$((found + 1))
    echo "[$found] Processing: $TREEFILE"

    python3 "$PYTHON_SCRIPT" "$TREEFILE" \
        --sh-threshold "$SH_THRESHOLD" \
        --bs-threshold "$BS_THRESHOLD"

    echo "-----------------------------------"

done < <(find "$INPUT_DIR" -type f -name "*.treefile" -print0)

if [ "$found" -eq 0 ]; then
    echo "Warning: no .treefile files found under '$INPUT_DIR'."
    exit 1
fi

echo "Done. Processed $found treefile(s)."