#!/usr/bin/env bash
# =============================================================================
# GraphAllele - example run script
#
#   bash run_GraphAllele.sh
#
# Edit the USER CONFIGURATION block below, then run from anywhere. The pipeline
# is breakpoint-resumable: re-running skips any chromosome-group step that has
# already completed successfully.
#
# HPC schedulers: this is a plain shell script with no scheduler directives.
# To submit it, wrap it in your scheduler's job script, e.g.:
#   PBS  -> prepend "#PBS -l nodes=1:ppn=12" lines, then `qsub run_GraphAllele.sh`
#   SLURM-> prepend "#SBATCH --cpus-per-task=12" lines, then `sbatch run_GraphAllele.sh`
# Keep the configured THREADS consistent with the cores you request.
# =============================================================================

set -eo pipefail

# Always run from the directory that contains this script (the workflow/ dir).
cd "$(dirname "$(readlink -f "$0")")"

# ======================= USER CONFIGURATION (edit me) ========================
# Conda environment created from environment.yml
ENV_NAME="polyalleler"

# --- Input files (relative to workflow/, or use absolute paths) --------------
GFF="../data/genome.gff3"            # whole-genome GFF3 annotation
FASTA="../data/genome.fasta"         # whole-genome FASTA sequence
REF_GFF="../data/reference.gff3"     # reference genome GFF3 (for calibration)
REF_CDS="../data/reference.cds"      # reference genome CDS FASTA (for TBLASTN)

# --- Orthogroups: choose ONE option -----------------------------------------
#   AUTO_OG=1 : run intra-group OrthoFinder automatically (ORTHOGROUPS ignored)
#   AUTO_OG=0 : supply a pre-computed Orthogroups.tsv via ORTHOGROUPS
AUTO_OG=1
ORTHOGROUPS=""                       # e.g. ../data/Orthogroups.tsv

# --- Chromosome range and subgenome layout ----------------------------------
START=1                              # first chromosome group number
END=10                               # last chromosome group number
# Haplotype suffixes used in your chromosome names Chr{NUM}{SUFFIX}.
# Must match exactly, e.g. for Chr1A..Chr1K use "A,B,C,D,E,F,G,H,I,J,K".
SUB_LIST="A,B,C,D,E,F,G,H,I,J,K,L,M,N"

# --- Numeric parameters ------------------------------------------------------
THREADS=12
MIN_C=3                              # min haplotypes a cluster must span
TANDEM_DIST=5                        # max gene-index distance for tandem detection
CLUSTER_DIST=30                      # max gene-index distance for synteny clustering
VERIFY_RATIO=0.35                    # min OrthoGroup support ratio for verification

# --- Output ------------------------------------------------------------------
OUTDIR="result/GraphAllele_out"
# =============================================================================


# --- activate conda environment (no hardcoded user paths) --------------------
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "[ERROR] conda not found. Install Miniconda/Anaconda or edit this script." >&2
    exit 1
fi
conda activate "$ENV_NAME"

# --- pre-flight: confirm required inputs exist and are non-empty -------------
echo "[INFO] Checking input files..."
REQUIRED=("$GFF" "$FASTA" "$REF_GFF" "$REF_CDS")
if [[ "$AUTO_OG" -ne 1 ]]; then
    REQUIRED+=("$ORTHOGROUPS")
fi
for f in "${REQUIRED[@]}"; do
    if [[ ! -s "$f" ]]; then
        echo "[ERROR] Required input missing or empty: $f" >&2
        exit 1
    fi
    echo "  [OK] $f  ->  $(readlink -f "$f")"
done

# --- assemble the argument list ----------------------------------------------
ARGS=(
    -g "$GFF"
    -f "$FASTA"
    -ref_g "$REF_GFF"
    -ref_f "$REF_CDS"
    -s "$START" -e "$END"
    -t "$THREADS"
    -o "$OUTDIR"
    --min_c "$MIN_C"
    --tandem_dist "$TANDEM_DIST"
    --cluster_dist "$CLUSTER_DIST"
    --verify_ratio "$VERIFY_RATIO"
    --sub_list "$SUB_LIST"
)
if [[ "$AUTO_OG" -eq 1 ]]; then
    ARGS+=(--auto_og)
else
    ARGS+=(-og "$ORTHOGROUPS")
fi

# --- run ---------------------------------------------------------------------
echo "[INFO] Launching GraphAllele (Chr${START}-${END}, ${THREADS} threads)..."
python GraphAllele.py "${ARGS[@]}"

echo "[INFO] GraphAllele finished. Results in: $OUTDIR"
