# GraphAllele

A graph-constrained pipeline for constructing standardized allele matrices in polyploid genomes.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Pipeline Architecture](#pipeline-architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Output Description](#output-description)
- [Dependencies](#dependencies)
- [Changelog](#changelog)
- [License](#license)
- [Contact](#contact)

---

## Overview

GraphAllele integrates synteny-based graph clustering, tandem duplication filtering, orthogroup verification, and reference genome calibration into a unified, resumable pipeline for allele identification in complex polyploid genomes.

**Chromosome naming convention required**: `Chr{NUM}{SUFFIX}`, e.g., `Chr1A`, `Chr1B`, ..., `Chr1K`.

---

## Features

- **Breakpoint-resumable**: each step is checkpointed and skipped automatically on re-run; if the pipeline fails at any step, simply re-run the same command and it will resume from the last incomplete step.
- **Sequential processing**: chromosome groups are processed one by one to protect HPC NFS I/O and ensure stable, readable logging.
- **Graph-constrained clustering**: NetworkX-based extraction of syntenic connected components with configurable distance thresholds.
- **Tandem duplication filtering**: self-BLASTP + graph-based tandem array blacklisting with adjustable gene distance.
- **Intra-group OrthoFinder with early termination**: OrthoFinder is run independently per chromosome group; the process is killed immediately after `Orthogroups.tsv` is written to disk, skipping the MSA and tree-building phases to reduce runtime.
- **Reference calibration**: TBLASTN-based anchoring to an external reference genome for cross-species annotation.
- **Standardized output**: fixed-column allele matrix (A–N subgenomes) with globally unique cluster IDs.
- **Automated JCVI ortholog analysis**: built-in shell script for pairwise synteny comparisons with auto-discovery of species from input files.
- **Global matrix normalization**: post-processing module that produces both a cleaned global matrix and a k-mer–formatted cluster file for downstream analysis.

---

## Pipeline Architecture

```
Input: Whole-genome GFF3 + FASTA + Reference CDS/GFF + Orthogroups.tsv (or --auto_og)
         |
         v
  Step 1: Prepare Data (01.prepare/)
          Split GFF/FASTA by chromosome group
          Extract CDS, PEP, BED per haplotype (gffread)
         |
         v
  Step 1.5 (optional): Intra-Group OrthoFinder (--auto_og)
          Run OrthoFinder per chromosome group
          Early termination after Orthogroups.tsv is written
         |
         v
  Step 2: Tandem Duplication Identification (02.tandem/)
          Self-BLASTP on merged PEP
          Graph-based tandem array detection (--tandem_dist)
         |
         v
  Step 3: JCVI Synteny Analysis (03.jcvi/)
          Pairwise LAST alignment via run_ortholog.sh
          MCScan anchor generation (.anchors files)
         |
         v
  Step 4: Graph Clustering (04.cluster.tsv)
          Build synteny graph from .anchors
          Filter tandem blacklist
          Extract connected components as candidate allele clusters
         |
         v
  Step 5: OrthoGroup Verification (05.verified.tsv)
          Validate clusters against OrthoFinder orthogroups
          Output verified and rejected tables separately
         |
         v
  Step 6: Sequence Homology Expansion (06.expanded/)
          BLASTP-based refinement of allele assignments
         |
         v
  Step 7: Reference Calibration (07.FINAL_ALLELE.tsv)
          TBLASTN anchoring to reference CDS
          Annotate each cluster with Ref_Gene and Ref_Locus
         |
         v
  Post-processing: Global Merge & Normalization
          Merge all chromosome groups into PolyAlleler_Global_Matrix.tsv
          Normalize into PolyAlleler_Global_Matrix_Cleaned.tsv (fixed A–N columns)
          Generate my_clusters.tsv (k-mer format)
```

---

## Installation

### Conda (Recommended)

```bash
git clone https://github.com/GengruiZhu/GraphAllele.git
cd GraphAllele

conda env create -f environment.yml
conda activate polyalleler
```

### Manual Installation

See [INSTALL.md](INSTALL.md) for step-by-step instructions.

---

## Quick Start

### Basic run with a pre-computed Orthogroups file

```bash
cd GraphAllele/workflow

python GraphAllele.py \
  -g /path/to/genome.gff3 \
  -f /path/to/genome.fasta \
  -ref_g /path/to/reference.gff3 \
  -ref_f /path/to/reference.cds \
  -og /path/to/Orthogroups.tsv \
  -s 1 -e 5 \
  -t 20 \
  --sub_list A,B,C,D,E
```

### Run with automatic OrthoFinder (no pre-computed file needed)

```bash
python GraphAllele.py \
  -g /path/to/genome.gff3 \
  -f /path/to/genome.fasta \
  -ref_g /path/to/reference.gff3 \
  -ref_f /path/to/reference.cds \
  --auto_og \
  -s 1 -e 10 \
  -t 30 \
  --min_c 2 \
  --cluster_dist 100 \
  --sub_list A,B,C,D,E,F,G,H,I,J,K,L,M,N
```

### PBS Job Submission

```bash
# Edit workflow/run_GraphAllele.sh to set your paths and resource requirements, then:
qsub workflow/run_GraphAllele.sh
```

**Example PBS script** (`workflow/run_GraphAllele.sh`):

```bash
#!/bin/bash
#PBS -N GraphAllele_ZG
#PBS -l nodes=1:ppn=30
#PBS -q comput
#PBS -l mem=
#PBS -o Allele.log
#PBS -j oe

cd $PBS_O_WORKDIR
date -R
source ~/miniconda3/bin/activate
source activate polyalleler

python GraphAllele.py \
  -g ../data/ZG.gff3 \
  -f ../data/ZG.fasta \
  -ref_g ../data/Eru.gff3 \
  -ref_f ../data/Eru.cds \
  --auto_og \
  -s 1 -e 10 \
  -t 30 \
  -o ZG_100 \
  --min_c 2 \
  --cluster_dist 100 \
  --sub_list A,B,C,D,E,F,G,H,I,J,K,L,M,N
```

---

## Usage

```
python GraphAllele.py [OPTIONS]
```

### Required Arguments

| Argument | Description |
|---|---|
| `-g / --gff` | Whole-genome GFF3 annotation file. Must use `Chr{NUM}{SUFFIX}` naming convention. |
| `-f / --fasta` | Whole-genome FASTA sequence file corresponding to the GFF3. |
| `-ref_g / --ref_gff` | Reference genome GFF3 annotation (used for calibration in Step 7). |
| `-ref_f / --ref_cds` | Reference genome CDS FASTA (used for TBLASTN in Step 7). |

### Optional Arguments

| Argument | Default | Description |
|---|---|---|
| `-og / --orthogroups` | `None` | Path to a pre-computed OrthoFinder `Orthogroups.tsv`. If provided, the pipeline uses this file for orthogroup verification (Step 5). Mutually exclusive in spirit with `--auto_og`, though both can be set (pre-computed file takes priority if the auto-generated one does not exist). |
| `--auto_og` | `False` | Run OrthoFinder automatically per chromosome group using only the haplotype PEP files for that group. The OrthoFinder process is terminated immediately after `Orthogroups.tsv` is written, skipping MSA and tree-building. Recommended when no pre-computed orthogroups are available. |
| `-s / --start` | `1` | Start chromosome group number. For example, `-s 1` means beginning from `Group_Chr01`. |
| `-e / --end` | `10` | End chromosome group number. For example, `-e 10` means processing through `Group_Chr10`. |
| `-t / --threads` | `10` | Total CPU threads available. Used by BLAST, OrthoFinder, and other parallelizable steps. |
| `-o / --outdir` | `standardized_results` | Output directory. All intermediate and final results are stored here. |
| `--sub_list` | `A,B,...,N` | Comma-separated haplotype suffix list. Defines the column order in the output allele matrix. Must match the suffix letters in your chromosome names. |
| `--min_c` | `3` | Minimum number of distinct haplotype chromosomes a cluster must span to be retained. Lower values (e.g., 2) are more permissive and suitable for genomes with high gene loss. |
| `--tandem_dist` | `5` | Maximum gene index distance for tandem duplicate detection (MCScanX convention). Two genes within this distance on the same chromosome with high sequence similarity are flagged as tandem duplicates. |
| `--cluster_dist` | `30` | Maximum gene index distance for synteny graph clustering. Increasing this value (e.g., to 100) allows more distant syntenic gene pairs to be grouped together, which can be useful for genomes with many structural rearrangements. |

### Parameter Guidance

- **`--sub_list`**: must exactly match the suffix letters used in your chromosome names. For example, if your chromosomes are named `Chr1A` through `Chr1K`, use `--sub_list A,B,C,D,E,F,G,H,I,J,K`. The order of suffixes determines column order in the output matrix.
- **`--min_c`**: a cluster must span at least this many distinct haplotype chromosomes to be retained. For highly fragmented or recently diverged polyploids, consider lowering to `2`.
- **`--auto_og`**: when enabled, OrthoFinder is run once per chromosome group using only the haplotype PEP files for that group. The process is terminated automatically after `Orthogroups.tsv` is written, before MSA and tree steps. It is still recommended to supply a pre-computed `-og` file when one is available, as it tends to be more comprehensive.
- **`--cluster_dist`**: the default value of `30` works well for most plant genomes. For genomes with extensive rearrangements (e.g., sugarcane), values of `50–100` may produce better results.
- **`--tandem_dist`**: the default of `5` follows the MCScanX convention. Increase to `10` if you suspect many dispersed tandem arrays are being missed.

---

## Output Description

```
standardized_results/
├── Group_Chr01/
│   ├── 01.prepare/                    # Per-haplotype GFF, FASTA, CDS, BED, PEP
│   ├── 01.5.OrthoFinder_Intra/        # Intra-group OrthoFinder results (--auto_og only)
│   ├── 02.tandem/                     # Merged GFF/PEP + tandem blacklist
│   ├── 03.jcvi/                       # Pairwise .anchors files
│   ├── 04.cluster.tsv                 # Raw synteny clusters
│   ├── 05.verified.tsv                # OrthoGroup-validated clusters
│   ├── 05.verified_rejected.tsv       # Clusters failing OG verification
│   ├── 06.expanded_expanded.tsv       # BLAST-refined allele table
│   └── 07.FINAL_ALLELE.tsv            # Final allele matrix with reference anchors
├── Group_Chr02/
│   └── ...
├── PolyAlleler_Global_Matrix.tsv       # Genome-wide merged allele matrix (raw)
├── PolyAlleler_Global_Matrix_Cleaned.tsv  # Normalized matrix with fixed A–N columns
└── my_clusters.tsv                     # K-mer format cluster file
```

### `07.FINAL_ALLELE.tsv` column format

| ClusterID | Ref_Gene | Ref_Locus | Chr01A | Chr01B | ... |
|---|---|---|---|---|---|
| Global_Cluster_000001 | Gene001 | Chr1:1000-5000(+) | gene_A | gene_B | ... |

Each row is one allele group. `-` or `NA` indicates no gene assigned for that haplotype in this cluster.

### `PolyAlleler_Global_Matrix_Cleaned.tsv`

A post-processed version of the global matrix with standardized columns (`Allele_A` through `Allele_N`) and globally re-numbered cluster IDs, produced by the normalization module.

### `my_clusters.tsv`

A two-column file: `ClusterID` and a comma-separated list of all gene members in that cluster. Useful for downstream k-mer or expression analyses.

---

## Dependencies

### Python Packages

| Package | Version |
|---|---|
| Python | >= 3.8 |
| Biopython | >= 1.79 |
| pandas | >= 1.3 |
| networkx | >= 2.6 |

### External Tools

| Tool | Purpose |
|---|---|
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | Sequence similarity search (BLASTP, TBLASTN, makeblastdb) |
| [gffread](https://github.com/gpertea/gffread) | CDS/protein extraction from GFF |
| [JCVI/MCScan](https://github.com/tanghaibao/jcvi) | Synteny analysis |
| [LAST](https://gitlab.com/mcfrith/last) | Pairwise sequence alignment |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) *(optional)* | Orthogroup inference (required if using `--auto_og`) |

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a detailed version history.

### v1.1 (2025-04-02)

- Added automated JCVI ortholog shell script (`bin/run_ortholog.sh`) with species auto-discovery
- Added global matrix normalization module (`generate_clean_clusters_auto`) producing cleaned matrix and k-mer cluster file
- Improved global merge logic using pandas with re-numbered ClusterIDs
- Added `--protein` flag to `prepare_jcvi.py` for protein sequence extraction
- Enhanced error handling with `[FATAL ERROR]` / `[ACTION]` messages and resume guidance
- Updated PBS example script with real-world parameters

### v1.0

- Initial release

---

## Citation

> Manuscript in preparation. Citation will be added upon publication.

---

## License

This project is licensed under a **Non-Commercial Research License**.
Free to use for academic and research purposes only. Commercial use is strictly prohibited.
See the [LICENSE](LICENSE) file for details.

---

## Contact

- **Developers**: Gengrui Zhu, Yi Chen
- **Issues**: please use the [GitHub Issues](https://github.com/GengruiZhu/GraphAllele/issues) page for bug reports and questions.
