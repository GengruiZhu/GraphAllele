#!/bin/bash
INPUT_DIR="./"
CSCORE=0.99

species=()
for cds_file in "$INPUT_DIR"/*.cds; do
    base=$(basename "$cds_file" .cds)
    bed_file="$INPUT_DIR/$base.bed"
    if [[ -f "$bed_file" ]]; then
        species+=("$base")
    fi
done

if [[ ${#species[@]} -lt 2 ]]; then
    echo "[ERROR] Less than 2 valid species found!"
    exit 1
fi

for (( i=0; i<${#species[@]}; i++ )); do
    for (( j=i+1; j<${#species[@]}; j++ )); do
        a="${species[$i]}"
        b="${species[$j]}"
        echo "[INFO] Running JCVI ortholog: $a vs $b"
        python -m jcvi.compara.catalog ortholog --no_strip_names "$a" "$b" --cscore=$CSCORE --no_dotplot
    done
done

echo "[DONE] JCVI ortholog analysis completed!"
