# Changelog

All notable changes to GraphAllele are documented here.

## v1.2

Maintenance and packaging release. The core algorithms are unchanged from the
internal development build; this release standardizes the package layout and
removes environment-specific code so it can be run anywhere.

### Packaging
- Consolidated the per-dataset run scripts into a single, fully parameterized `workflow/run_GraphAllele.sh` with a `USER CONFIGURATION` block at the top.
- Removed all hardcoded absolute paths, user accounts, and references to specific conda installations. Conda activation now resolves the environment via `conda info --base` with `$HOME`-based fallbacks.
- Removed scheduler (PBS) directives from the run script; it is now a plain shell script that can be wrapped for any scheduler.
- Removed internal/experimental scripts that were not part of the pipeline (batch parameter-sweep submitter, per-dataset result collectors, and standalone evaluation scripts).
- Added `data/README.md` describing the expected inputs.
- Updated `README.md`, `INSTALL.md`, `environment.yml`, and `.gitignore` to match the actual code and outputs.

### Documentation fixes
- Documented the `--verify_ratio` parameter (default 0.35).
- Corrected the documented output filenames: the genome-wide matrix is `PolyAlleler_Global_Matrix_Cleaned_<paramtag>.tsv`, accompanied by a flat `my_clusters_<paramtag>.tsv` list.
- Corrected the documented column layout: the genome-wide matrix uses `Allele_<SUFFIX>` columns plus a `Salvaged_Genes` column, while the per-group `07.FINAL_ALLELE.tsv` uses `Chr{NN}{SUFFIX}` columns plus `Unplaced_or_Translocated`.

### Algorithms (carried over from the internal build)
- **Unified gene-ID handling** (`safe_extract_id`): consistent stripping of `gene:` / `transcript:` prefixes and trailing transcript suffixes across all modules.
- **BLAST sandboxing**: all BLAST databases are built inside temporary directories to avoid polluting reference data directories.
- **Functional sequence-homology expansion**: BLASTP now recovers high-identity unclustered paralogs/translocated copies into the matrix.
- **`Unplaced_or_Translocated` handling** in clustering and verification, with comma-separated multi-gene cells.
- **OrthoGroup-based gene rescue** in verification, placing same-OG members back into the correct subgenome column.
- **Tactical Sniper** early termination for intra-group OrthoFinder: backs up `Orthogroups.tsv` on detection, then SIGTERM → (grace period) → SIGKILL of the process group, skipping MSA/tree steps.
- **Smart Subgenome Backfilling** when merging the genome-wide matrix.
- **Breakpoint resume** for tandem and synteny steps.

## v1.0

- Initial release.
