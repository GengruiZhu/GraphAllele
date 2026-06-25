# Installation

GraphAllele is a Python pipeline that orchestrates several external
bioinformatics tools. The recommended way to install everything is via Conda.

---

## 1. Conda (recommended)

```bash
git clone https://github.com/GengruiZhu/GraphAllele.git
cd GraphAllele

conda env create -f environment.yml
conda activate polyalleler
```

This creates an environment named `polyalleler` containing Python and all
required packages and external tools (BLAST+, gffread, JCVI, LAST, and
OrthoFinder).

Verify the key tools are visible:

```bash
python   --version
makeblastdb -version
blastp   -version
tblastn  -version
gffread  --version
python -m jcvi.compara.catalog --help   >/dev/null && echo "jcvi OK"
orthofinder -h | head -n 1               # only needed for --auto_og
```

---

## 2. Manual installation

If you prefer not to use the provided environment file, install the
dependencies yourself.

### 2.1 Python packages

Python >= 3.8 with:

```bash
pip install "biopython>=1.79" "pandas>=1.3" "networkx>=2.6"
```

### 2.2 External tools

Install these and make sure they are on your `PATH`:

| Tool | Notes |
| --- | --- |
| BLAST+ (>= 2.12) | provides `makeblastdb`, `blastp`, `tblastn` |
| gffread | CDS / protein extraction from GFF3 |
| LAST | pairwise aligner used by JCVI/MCScan |
| JCVI | `pip install jcvi` (requires LAST on PATH) |
| OrthoFinder | optional; only needed when running with `--auto_og` |

The fastest route for the external tools is still Bioconda:

```bash
conda install -c conda-forge -c bioconda blast last gffread jcvi orthofinder
```

---

## 3. Test the installation

```bash
cd GraphAllele/workflow
python GraphAllele.py --help
```

If the help text prints without import errors, the Python side is ready. Then
edit `run_GraphAllele.sh` (or pass arguments directly) to point at your data
and launch a run. See the [README](README.md) for usage details.

---

## 4. Notes

- The pipeline writes BLAST databases into temporary sandboxes, so it does not
  pollute your reference data directory.
- Each step is checkpointed; re-running the same command resumes from where it
  stopped rather than restarting from scratch.
- `run_GraphAllele.sh` activates the `polyalleler` environment automatically via
  `conda info --base` (falling back to `$HOME/miniconda3` or `$HOME/anaconda3`),
  so it contains no hardcoded, account-specific paths.
