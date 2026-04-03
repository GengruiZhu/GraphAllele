\# Changelog



All notable changes to GraphAllele are documented in this file.



\---



\## \[v1.1] - 2025-04-02



\### Added



\- \*\*`bin/run\_ortholog.sh`\*\*: New automated JCVI pairwise ortholog analysis script.

&#x20; - Automatically discovers species from `.cds` and `.bed` file pairs in the input directory.

&#x20; - Iterates over all pairwise species combinations and runs `jcvi.compara.catalog ortholog`.

&#x20; - Configurable `CSCORE` threshold (default: 0.99) and `--no\_strip\_names` mode.

&#x20; - Validates that at least 2 valid species are found before proceeding.



\- \*\*`generate\_clean\_clusters\_auto()` function\*\* in the main workflow:

&#x20; - Post-processing module that merges all `Group\_Chr\*/07.FINAL\_ALLELE.tsv` files into a single normalized global matrix.

&#x20; - Outputs `PolyAlleler\_Global\_Matrix\_Cleaned.tsv` with standardized columns (`ClusterID`, `Ref\_Gene`, `Ref\_Locus`, `Allele\_A` through `Allele\_N`).

&#x20; - Outputs `my\_clusters.tsv` in k-mer format (`ClusterID` + comma-separated gene list) for downstream analysis.

&#x20; - Handles multi-value allele entries (comma-separated genes within a cell) and maps them to the correct subgenome columns.



\- \*\*`--protein` flag\*\* added to the `prepare\_jcvi.py` call in Step 1 to ensure protein sequences are extracted alongside CDS.



\### Changed



\- \*\*Global merge logic improved\*\*: The pipeline now uses `pandas.concat()` to merge all per-group `07.FINAL\_ALLELE.tsv` files, with globally re-numbered `ClusterID` values (`Global\_Cluster\_000000`, `Global\_Cluster\_000001`, ...). After the raw merge, the new `generate\_clean\_clusters\_auto()` function is called to produce the normalized output.



\- \*\*Sequential execution messaging\*\*: Pipeline startup now prints `"Starting pipeline sequentially to protect HPC NFS I/O..."` to explicitly communicate that chromosome groups are processed one by one, which avoids I/O contention on HPC shared filesystems.



\- \*\*Error handling enhanced\*\*:

&#x20; - On failure, the pipeline now prints `\[FATAL ERROR]` with the specific chromosome group that failed.

&#x20; - An `\[ACTION]` message guides the user to check error logs, fix the issue, and re-run (the pipeline will automatically resume from the failed group).

&#x20; - The pipeline exits with `sys.exit(1)` on failure instead of silently continuing.



\- \*\*PBS example script\*\* (`workflow/run\_GraphAllele.sh`) updated with a real-world example using `--auto\_og`, `--cluster\_dist 100`, `--min\_c 2`, and a full A–N subgenome list.



\### Fixed



\- Single-group edge case: when only one chromosome group is processed (e.g., `-s 5 -e 5`), the pipeline now correctly skips the merge step and reports `"Just Single Group, skipping merging."`.



\---



\## \[v1.0] - 2025-03-20



\### Added



\- Initial release of GraphAllele.

\- 7-step pipeline: Prepare → Tandem → JCVI Synteny → Graph Clustering → OrthoGroup Verification → BLAST Expansion → Reference Calibration.

\- Breakpoint-resumable execution with per-step checkpointing.

\- Sequential chromosome group processing.

\- Intra-group OrthoFinder with early termination (Visual Sniper).

\- Support for up to 14 subgenomes (A–N).

\- PBS job submission template.

\- Conda environment file (`environment.yml`).

