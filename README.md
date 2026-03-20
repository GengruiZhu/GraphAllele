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
- [License](#license)
- [Contact](#contact)

---

## Overview

GraphAllele integrates synteny-based graph clustering, tandem duplication filtering, orthogroup verification, and reference genome calibration into a unified, resumable pipeline for allele identification in complex polyploid genomes.

Chromosome naming convention required: `Chr{NUM}{SUFFIX}`, e.g., `Chr1A`, `Chr1B`, ..., `Chr1K`.

---

## Features

- **Breakpoint-resumable**: each step is checkpointed and skipped automatically on re-run
- **Parallel processing**: chromosome groups are processed concurrently via Python multiprocessing
- **Graph-constrained clustering**: NetworkX-based extraction of syntenic connected components
- **Tandem duplication filtering**: self-BLASTP + graph-based tandem array blacklisting
- **OrthoFinder integration**: supports pre-computed `Orthogroups.tsv` or auto-run mode
- **Reference calibration**: TBLASTN-based anchoring to an external reference genome
- **Standardized output**: fixed-column allele matrix with globally unique cluster IDs

---

## Pipeline Architecture

```
Input: Whole-genome GFF3 + FASTA + Reference CDS/GFF + Orthogroups.tsv
         |
         v
  Step 1: Prepare Data
          Split GFF/FASTA by chromosome group
          Extract CDS, PEP, BED per haplotype (gffread)
         |
         v
  Step 2: Tandem Duplication Identification
          Self-BLASTP on merged PEP
          Graph-based tandem array detection (max_dist, min_identity)
         |
         v
  Step 3: JCVI Synteny Analysis
          Pairwise LAST alignment
          MCScan anchor generation (.anchors files)
         |
         v
  Step 4: Graph Clustering
          Build synteny graph from .anchors
          Filter tandem blacklist
          Extract connected components as candidate allele clusters
         |
         v
  Step 5: OrthoGroup Verification
          Validate clusters against OrthoFinder orthogroups
          Output verified and rejected tables separately
         |
         v
  Step 6: Sequence Homology Expansion
          BLASTP-based refinement of allele assignments
         |
         v
  Step 7: Reference Calibration
          TBLASTN anchoring to reference CDS
          Annotate each cluster with Ref_Gene and Ref_Locus
         |
         v
Output: 07.FINAL_ALLELE.tsv
        PolyAlleler_Global_Matrix.tsv  (multi-chromosome merged)
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
# Edit workflow/run_GraphAllele.sh to set your paths and resource requirements, then:
qsub workflow/run_GraphAllele.sh
```

---

## Usage

```
python GraphAllele.py [OPTIONS]
```

### Required Arguments

| Argument | Description |
|---|---|
| `-g / --gff` | Whole-genome GFF3 annotation file |
| `-f / --fasta` | Whole-genome FASTA sequence file |
| `-ref_g / --ref_gff` | Reference genome GFF3 (for calibration) |
| `-ref_f / --ref_cds` | Reference genome CDS FASTA (for TBLASTN) |

### Optional Arguments

| Argument | Default | Description |
|---|---|---|
| `-og / --orthogroups` | None | Pre-computed OrthoFinder `Orthogroups.tsv` |
| `--auto_og` | False | Auto-run OrthoFinder internally (slow for large genomes) |
| `-s / --start` | 1 | Start chromosome number |
| `-e / --end` | 10 | End chromosome number |
| `-t / --threads` | 10 | Total CPU threads |
| `-o / --outdir` | `standardized_results` | Output directory |
| `--sub_list` | `A,B,...,N` | Comma-separated haplotype suffix list |
| `--min_c` | 3 | Minimum haplotypes required per allele cluster |
| `--tandem_dist` | 5 | Max gene index distance for tandem detection |
| `--cluster_dist` | 30 | Max gene index distance for synteny clustering |

### Notes

**`--sub_list`**: must match the suffix letters used in your chromosome names. For example, if your chromosomes are named `Chr1A` through `Chr1K`, use `--sub_list A,B,C,D,E,F,G,H,I,J,K`.

**`--min_c`**: a cluster must span at least this many distinct haplotype chromosomes to be retained.

**OrthoFinder**: it is strongly recommended to run OrthoFinder independently and supply the result via `-og`. The `--auto_og` mode is provided for convenience but is not suitable for large polyploid genomes.

---

## Output Description

```
standardized_results/
├── Group_Chr01/
│   ├── 01.prepare/               # Per-haplotype GFF, FASTA, CDS, BED, PEP
│   ├── 02.tandem/                # Merged GFF/PEP + tandem blacklist
│   ├── 03.jcvi/                  # Pairwise .anchors files
│   ├── 04.cluster.tsv            # Raw synteny clusters
│   ├── 05.verified.tsv           # OrthoGroup-validated clusters
│   ├── 05.verified_rejected.tsv  # Clusters failing OG verification
│   ├── 06.expanded_expanded.tsv  # BLAST-refined allele table
│   └── 07.FINAL_ALLELE.tsv       # Final allele matrix with reference anchors
├── Group_Chr02/
│   └── ...
└── PolyAlleler_Global_Matrix.tsv  # Genome-wide merged allele matrix
```

### `07.FINAL_ALLELE.tsv` column format

| ClusterID | Ref_Gene | Ref_Locus | Chr01A | Chr01B | ... |
|---|---|---|---|---|---|
| Global_Cluster_000001 | Gene001 | Chr1:1000-5000(+) | gene_A | gene_B | ... |

Each row is one allele group. `NA` indicates no gene assigned for that haplotype in this cluster.

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
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | Sequence similarity search |
| [gffread](https://github.com/gpertea/gffread) | CDS/protein extraction from GFF |
| [JCVI/MCScan](https://github.com/tanghaibao/jcvi) | Synteny analysis |
| [LAST](https://gitlab.com/mcfrith/last) | Pairwise sequence alignment |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) *(optional)* | Orthogroup inference |

---
## Citation

> Manuscript in preparation. Citation will be added upon publication.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

- **Developers**: Gengrui Zhu, Yi Chen
- **Issues**: please use the [GitHub Issues](https://github.com/GengruiZhu/GraphAllele/issues) page for bug reports and questions.
