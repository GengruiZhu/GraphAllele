#!/bin/bash
# Enable strict error handling: script will exit if any unhandled command fails
set -e

INPUT_DIR="./"
CSCORE=0.99

# CRITICAL FIX: Prevent literal string expansion if no .cds files are found
shopt -s nullglob

species=()
# Step 1: Scan and validate input pairs
for cds_file in "$INPUT_DIR"/*.cds; do
    base=$(basename "$cds_file" .cds)
    bed_file="$INPUT_DIR/$base.bed"

    # Ensure both sequence and coordinate files exist
    if [[ -f "$bed_file" ]]; then
        species+=("$base")
    else
        echo "[WARNING] Missing corresponding BED file for $base, skipping this subgenome."
    fi
done

# Step 2: Ensure enough valid datasets
if [[ ${#species[@]} -lt 2 ]]; then
    echo "[ERROR] Less than 2 valid species/subgenomes found in $INPUT_DIR!"
    echo "[ERROR] Cannot perform pairwise synteny analysis. Exiting."
    exit 1
fi

echo "[INFO] Successfully loaded ${#species[@]} valid subgenomes. Initiating combinatorial pairwise alignments..."

# Step 3: Perform all-vs-all pairwise alignments
for (( i=0; i<${#species[@]}; i++ )); do
    for (( j=i+1; j<${#species[@]}; j++ )); do
        a="${species[$i]}"
        b="${species[$j]}"

        # TARGET ANCHOR CHECK (Breakpoint Resume Mechanism)
        # JCVI generates an output file named a.b.anchors. We check if it exists and is not empty.
        target_anchors="${a}.${b}.anchors"
        reverse_anchors="${b}.${a}.anchors" # Just in case JCVI sorted the names differently

        if [[ -s "$target_anchors" ]] || [[ -s "$reverse_anchors" ]]; then
            echo "[SKIP] $a vs $b: Anchors file already exists. Skipping computation."
            continue
        fi

        echo "[INFO] Running JCVI synteny pipeline: $a vs $b (C-score threshold: $CSCORE)"

        # Disable strict exit temporarily to catch the JCVI specific error and log it properly
        set +e
        python -m jcvi.compara.catalog ortholog --no_strip_names "$a" "$b" --cscore=$CSCORE --no_dotplot
        exit_code=$?
        set -e # Re-enable strict exit

        # Step 4: Halt the entire pipeline if a pairwise alignment crashes
        if [[ $exit_code -ne 0 ]]; then
            echo "[FATAL ERROR] JCVI crashed during the alignment of $a vs $b!"
            echo "[ACTION] Please check the LAST output logs above for out-of-memory (OOM) or format errors."
            exit 1
        fi
    done
done

echo "[SUCCESS] All pairwise JCVI ortholog analyses have been completed safely!"
