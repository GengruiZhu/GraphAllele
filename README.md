# GraphAllele

> **Graph-Constrained Allele Matrix Pipeline for Polyploid Genomes**
>
> A high-performance, modular bioinformatics pipeline for identifying and constructing standardized allele matrices in complex polyploid genomes (e.g., sugarcane, wheat, cotton).

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Pipeline Architecture](#pipeline-architecture)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Parameters](#parameters)
- [Output Description](#output-description)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [License](#license)

---

## Overview

GraphAllele (also known as **PolyAlleler V6.2**) is designed to handle the complexity of high-ploidy genomes by integrating synteny-based graph clustering, tandem duplication filtering, orthogroup verification, and reference genome calibration into a unified, resumable pipeline.

The pipeline was developed by engineers from **Zhangqing-Lab** (Yi Chen, Gengrui Zhu, et al.) and is particularly suited for highly polyploid species such as *Saccharum officinarum* (sugarcane) and related Poaceae.

---

## Features

- ✅ **Breakpoint-resumable**: Each step is checkpointed; interrupted runs automatically resume from the last completed step
- ✅ **Parallel processing**: Chromosome groups are processed concurrently via Python multiprocessing
- ✅ **Graph-constrained clustering**: Uses NetworkX graph algorithms to identify syntenic allele groups
- ✅ **Tandem duplication filtering**: Rigorous BLAST + graph-based identification of tandem arrays
- ✅ **OrthoFinder integration**: Supports both pre-computed and auto-generated orthogroups
- ✅ **Reference calibration**: TBLASTN-based anchoring to an external reference genome
- ✅ **Standardized output**: Fixed-column allele matrix with global unique cluster IDs

---

## Pipeline Architecture

```
Input: Whole-genome GFF3 + FASTA + Reference CDS/GFF + Orthogroups.tsv
         │
         ▼
  Step 1: Prepare Data
  (Split by chromosome group, extract CDS/PEP/BED per haplotype)
         │
         ▼
  Step 2: Tandem Duplication Identification
  (Self-BLASTP + Graph-based tandem array detection)
         │
         ▼
  Step 3: JCVI Synteny Analysis
  (Pairwise LAST alignment + MCScan anchor generation)
         │
         ▼
  Step 4: Graph Clustering
  (Build synteny graph, filter tandems, extract connected components)
         │
         ▼
  Step 5: OrthoGroup Verification
  (Validate clusters against OrthoFinder orthogroups)
         │
         ▼
  Step 6: Sequence Homology Expansion
  (BLASTP-based refinement of allele assignments)
         │
         ▼
  Step 7: Reference Calibration
  (TBLASTN anchoring to reference genome, locus annotation)
         │
         ▼
Output: 07.FINAL_ALLELE.tsv (Standardized Allele Matrix)
        PolyAlleler_Global_Matrix.tsv (Multi-chromosome merged matrix)
```

---

## Installation

### Option 1: Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/GraphAllele.git
cd GraphAllele

# Create and activate the conda environment
conda env create -f environment.yml
conda activate polyalleler
```

### Option 2: Manual Installation

See [INSTALL.md](INSTALL.md) for step-by-step manual installation instructions.

---

## Quick Start

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

### PBS Job Submission

```bash
# Edit workflow/run_GraphAllele.sh with your paths, then:
qsub workflow/run_GraphAllele.sh
```

---

## Usage

### Main Script

```
python GraphAllele.py [OPTIONS]
```

### Required Arguments

| Argument | Description |
|---|---|
| `-g / --gff` | Path to whole-genome GFF3 annotation file |
| `-f / --fasta` | Path to whole-genome FASTA sequence file |
| `-ref_g / --ref_gff` | Path to reference genome GFF3 (for calibration) |
| `-ref_f / --ref_cds` | Path to reference genome CDS FASTA (for TBLASTN) |

### Optional Arguments

| Argument | Default | Description |
|---|---|---|
| `-og / --orthogroups` | None | Pre-computed OrthoFinder `Orthogroups.tsv` |
| `--auto_og` | False | Auto-run OrthoFinder (not recommended for large genomes) |
| `-s / --start` | 1 | Start chromosome number |
| `-e / --end` | 10 | End chromosome number |
| `-t / --threads` | 10 | Number of CPU threads |
| `-o / --outdir` | `standardized_results` | Output directory |
| `--sub_list` | `A,B,...,N` | Comma-separated haplotype suffix list |
| `--min_c` | 3 | Minimum number of haplotypes per allele cluster |
| `--tandem_dist` | 5 | Maximum gene index distance for tandem detection |
| `--cluster_dist` | 30 | Maximum gene index distance for synteny clustering |

---

## Parameters

### `--sub_list`
Defines the expected haplotype suffixes for your polyploid genome. For example:
- Sugarcane (11-ploid A–K): `--sub_list A,B,C,D,E,F,G,H,I,J,K`
- Wheat (hexaploid A–C subgenomes): `--sub_list A,B,D`

Chromosomes must be named in the format `Chr{NUM}{SUFFIX}` in your GFF3/FASTA, e.g., `Chr1A`, `Chr1B`, ..., `Chr1K`.

### `--min_c`
Controls the strictness of allele cluster filtering. A cluster must span at least `--min_c` different haplotype chromosomes to be retained.

### OrthoFinder Integration
It is **strongly recommended** to run OrthoFinder independently before running GraphAllele, especially for large polyploid genomes where auto-OrthoFinder may take days to weeks. Provide the result via `-og Orthogroups.tsv`.

---

## Output Description

```
standardized_results/
├── Group_Chr01/
│   ├── 01.prepare/          # Per-haplotype GFF, FASTA, CDS, BED, PEP
│   ├── 02.tandem/           # Tandem duplication results + merged GFF/PEP
│   ├── 03.jcvi/             # Pairwise synteny anchors (.anchors files)
│   ├── 04.cluster.tsv       # Raw synteny clusters
│   ├── 05.verified.tsv      # OrthoGroup-validated clusters
│   ├── 05.verified_rejected.tsv  # Rejected clusters
│   ├── 06.expanded_expanded.tsv  # BLAST-refined allele table
│   └── 07.FINAL_ALLELE.tsv  # Final allele matrix with reference anchors
├── Group_Chr02/
│   └── ...
└── PolyAlleler_Global_Matrix.tsv  # Genome-wide merged allele matrix
```

### `07.FINAL_ALLELE.tsv` Format

| ClusterID | Ref_Gene | Ref_Locus | Chr01A | Chr01B | Chr01C | ... |
|---|---|---|---|---|---|---|
| Global_Cluster_000001 | EruGene001 | Chr1:1000-5000(+) | ROC_gene001 | ROC_gene002 | ROC_gene003 | ... |

---

## Dependencies

### Python Packages

| Package | Version |
|---|---|
| Python | ≥ 3.8 |
| Biopython | ≥ 1.79 |
| pandas | ≥ 1.3 |
| networkx | ≥ 2.6 |

### External Tools

| Tool | Purpose | Installation |
|---|---|---|
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | Sequence similarity search | `conda install -c bioconda blast` |
| [gffread](https://github.com/gpertea/gffread) | CDS/protein extraction from GFF | `conda install -c bioconda gffread` |
| [JCVI/MCScan](https://github.com/tanghaibao/jcvi) | Synteny analysis | `conda install -c bioconda jcvi` |
| [LAST](https://gitlab.com/mcfrith/last) | Pairwise sequence alignment | `conda install -c bioconda last` |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) *(optional)* | Orthogroup inference | `conda install -c bioconda orthofinder` |

---

## Citation

If you use GraphAllele in your research, please cite:

> Zhu G., Chen Y., et al. (2025). GraphAllele: A Graph-Constrained Pipeline for Allele Matrix Construction in Complex Polyploid Genomes. *Manuscript in preparation.*

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

- **Developer**: Gengrui Zhu, Yi Chen
- **Lab**: Zhangqing Lab
- **Issues**: Please open a GitHub Issue for bug reports and feature requests.
