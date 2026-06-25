# GraphAllele

A graph-constrained pipeline for constructing standardized allele matrices in
polyploid genomes.

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

GraphAllele integrates synteny-based graph clustering, tandem duplication
filtering, orthogroup verification, sequence-homology expansion, and reference
genome calibration into a unified, resumable pipeline for allele identification
in complex polyploid genomes.

Chromosome naming convention required: `Chr{NUM}{SUFFIX}`, e.g. `Chr1A`,
`Chr1B`, ..., `Chr1K`.

---

## Features

- **Breakpoint-resumable**: each step is checkpointed and skipped automatically on re-run.
- **Sequential processing**: chromosome groups are processed one by one for stable, readable logging and to protect HPC NFS I/O.
- **Graph-constrained clustering**: NetworkX-based extraction of syntenic connected components, with an intra-chromosomal gene-distance constraint.
- **Tandem duplication filtering**: self-BLASTP plus graph-based tandem-array blacklisting.
- **Intra-group OrthoFinder with early termination (Tactical Sniper)**: OrthoFinder is run independently per chromosome group; once `Orthogroups.tsv` is confirmed on disk, the result is backed up and the process group is shut down (SIGTERM, then SIGKILL after a short grace period), skipping the MSA and tree-building phases to save time.
- **Sequence-homology expansion**: BLASTP recovers high-identity unclustered paralogs/translocated copies into the matrix.
- **Smart Subgenome Backfilling**: salvaged genes whose IDs carry a subgenome signature are routed back into the correct allele column; only truly unplaceable genes land in the catch-all column.
- **Reference calibration**: TBLASTN-based anchoring to an external reference genome, annotating each cluster with `Ref_Gene` and `Ref_Locus`.
- **Standardized output**: fixed-column allele matrix with globally unique cluster IDs.

---

## Pipeline Architecture

```
Input: whole-genome GFF3 + FASTA + reference CDS/GFF
       (+ Orthogroups.tsv, or --auto_og)
         |
         v
  Step 1: Prepare Data
          Split GFF/FASTA by chromosome group
          Extract CDS, PEP, BED per haplotype (gffread)
         |
         v
  Step 1.5: Intra-Group OrthoFinder        (only with --auto_og)
          Per-group OrthoFinder + Tactical Sniper early termination
         |
         v
  Step 2: Tandem Duplication Identification
          Self-BLASTP on merged PEP
          Graph-based tandem-array detection (max_dist, min_identity)
         |
         v
  Step 3: JCVI Synteny Analysis
          Pairwise LAST alignment + MCScan anchors (.anchors)
         |
         v
  Step 4: Graph Clustering
          Build synteny graph from .anchors, drop tandem blacklist,
          extract connected components as candidate allele clusters
         |
         v
  Step 5: OrthoGroup Verification
          Validate clusters against OrthoFinder orthogroups,
          rescue same-OG members; write verified + rejected tables
         |
         v
  Step 6: Sequence Homology Expansion
          BLASTP recovery of high-identity unclustered genes
         |
         v
  Step 7: Reference Calibration
          TBLASTN anchoring to reference CDS; annotate Ref_Gene / Ref_Locus
         |
         v
Output: per-group 07.FINAL_ALLELE.tsv
        PolyAlleler_Global_Matrix_Cleaned_<paramtag>.tsv   (genome-wide merged)
        my_clusters_<paramtag>.tsv                          (flat gene list)
```

---

## Installation

### Conda (recommended)

```bash
git clone https://github.com/GengruiZhu/GraphAllele.git
cd GraphAllele

conda env create -f environment.yml
conda activate polyalleler
```

### Manual installation

See [INSTALL.md](INSTALL.md) for step-by-step instructions.

---

## Quick Start

Edit the `USER CONFIGURATION` block at the top of `workflow/run_GraphAllele.sh`
to point at your data, then run it directly:

```bash
cd GraphAllele/workflow
bash run_GraphAllele.sh
```

Or call the pipeline directly:

```bash
cd GraphAllele/workflow

python GraphAllele.py \
  -g ../data/genome.gff3 \
  -f ../data/genome.fasta \
  -ref_g ../data/reference.gff3 \
  -ref_f ../data/reference.cds \
  --auto_og \
  -s 1 -e 10 \
  -t 12 \
  --sub_list A,B,C,D,E,F,G,H,I,J,K \
  --min_c 3 \
  --cluster_dist 30 \
  --verify_ratio 0.35 \
  -o result/GraphAllele_out
```

### Running on an HPC scheduler

`run_GraphAllele.sh` is a plain shell script with no scheduler directives. To
submit it, wrap it in your scheduler's job script — e.g. prepend `#PBS` /
`#SBATCH` resource lines and submit with `qsub` / `sbatch`. Keep the configured
`THREADS` consistent with the cores you request.

---

## Usage

```bash
python GraphAllele.py [OPTIONS]
```

### Required arguments

| Argument | Description |
| --- | --- |
| `-g / --gff` | Whole-genome GFF3 annotation file |
| `-f / --fasta` | Whole-genome FASTA sequence file |
| `-ref_g / --ref_gff` | Reference genome GFF3 (for calibration) |
| `-ref_f / --ref_cds` | Reference genome CDS FASTA (for TBLASTN) |

### Orthogroups (supply one)

| Argument | Default | Description |
| --- | --- | --- |
| `-og / --orthogroups` | None | Pre-computed OrthoFinder `Orthogroups.tsv` |
| `--auto_og` | off | Run OrthoFinder automatically per chromosome group |

One of `-og` or `--auto_og` is required.

### Optional arguments

| Argument | Default | Description |
| --- | --- | --- |
| `-s / --start` | 1 | Start chromosome group number |
| `-e / --end` | 10 | End chromosome group number |
| `-t / --threads` | 10 | CPU threads |
| `-o / --outdir` | `standardized_results` | Output directory |
| `--sub_list` | `A,B,...,N` | Comma-separated haplotype suffix list |
| `--min_c` | 3 | Minimum haplotypes required per allele cluster |
| `--tandem_dist` | 5 | Max gene-index distance for tandem detection |
| `--cluster_dist` | 30 | Max gene-index distance for synteny clustering |
| `--verify_ratio` | 0.35 | Minimum OrthoGroup support ratio for cluster verification |

### Notes

**`--sub_list`** must match the suffix letters used in your chromosome names.
For example, if your chromosomes are named `Chr1A` through `Chr1K`, use
`--sub_list A,B,C,D,E,F,G,H,I,J,K`.

**`--min_c`**: a cluster must span at least this many distinct haplotype
chromosomes to be retained.

**`--verify_ratio`**: within a candidate cluster, the fraction of genes sharing
the dominant orthogroup must reach this value for the cluster to be verified.
Lower values are more permissive.

**`--auto_og`**: when enabled, OrthoFinder runs once per chromosome group using
only that group's haplotype PEP files, and is terminated automatically once
`Orthogroups.tsv` is written (before MSA and tree steps). Supplying a
pre-computed `-og` file is still recommended when one is available.

---

## Output Description

```
<outdir>/
├── Group_Chr01/
│   ├── 01.prepare/               # Per-haplotype GFF, FASTA, CDS, BED, PEP
│   ├── 01.5.OrthoFinder_Intra/   # Intra-group OrthoFinder results (--auto_og only)
│   ├── 02.tandem/                # Merged GFF/PEP + tandem blacklist
│   ├── 03.jcvi/                  # Pairwise .anchors files
│   ├── 04.cluster.tsv            # Raw synteny clusters
│   ├── 05.verified.tsv           # OrthoGroup-validated clusters
│   ├── 05.verified_rejected.tsv  # Clusters failing OG verification
│   ├── 06.expanded_expanded.tsv  # BLAST-refined allele table
│   └── 07.FINAL_ALLELE.tsv       # Per-group allele matrix with reference anchors
├── Group_Chr02/
│   └── ...
├── PolyAlleler_Global_Matrix_Cleaned_<paramtag>.tsv   # Genome-wide merged matrix
└── my_clusters_<paramtag>.tsv                          # Flat per-cluster gene list
```

`<paramtag>` encodes the run parameters, e.g. `minc3_dist30_ratio0.35`.

### Per-group `07.FINAL_ALLELE.tsv`

Columns: `ClusterID`, `Ref_Gene`, `Ref_Locus`, then one column per haplotype
named `Chr{NN}{SUFFIX}` (e.g. `Chr01A`, `Chr01B`, ...), plus
`Unplaced_or_Translocated`.

| ClusterID | Ref_Gene | Ref_Locus | Chr01A | Chr01B | ... | Unplaced_or_Translocated |
| --- | --- | --- | --- | --- | --- | --- |
| Cluster_01_00001 | Gene001 | Chr1:1000-5000(+) | gene_A | gene_B | ... | NA |

### Genome-wide `PolyAlleler_Global_Matrix_Cleaned_<paramtag>.tsv`

Columns: `ClusterID`, `Ref_Gene`, `Ref_Locus`, one `Allele_<SUFFIX>` column per
subgenome (e.g. `Allele_A`, `Allele_B`, ...), plus a final `Salvaged_Genes`
catch-all column for genes that could not be placed into a subgenome.

| ClusterID | Ref_Gene | Ref_Locus | Allele_A | Allele_B | ... | Salvaged_Genes |
| --- | --- | --- | --- | --- | --- | --- |
| Global_Cluster_000001 | Gene001 | Chr1:1000-5000(+) | gene_A | gene_B | ... | NA |

`my_clusters_<paramtag>.tsv` is a two-column flat list: `ClusterID` and a
comma-separated list of all genes in that cluster. `NA` indicates no gene
assigned for that haplotype/subgenome in a cluster.

---

## Dependencies

### Python packages

| Package | Version |
| --- | --- |
| Python | >= 3.8 |
| Biopython | >= 1.79 |
| pandas | >= 1.3 |
| networkx | >= 2.6 |

### External tools

| Tool | Purpose |
| --- | --- |
| BLAST+ | Sequence similarity search (makeblastdb, blastp, tblastn) |
| gffread | CDS/protein extraction from GFF |
| JCVI / MCScan | Synteny analysis |
| LAST | Pairwise sequence alignment (used by JCVI) |
| OrthoFinder *(optional)* | Orthogroup inference, required only for `--auto_og` |

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for the full history.

---

## License

This project is licensed under a **Non-Commercial Research License**. It is free
to use for academic and non-profit research purposes only; commercial use
requires a separate written agreement. See the [LICENSE](LICENSE) file for the
full text.

---

## Contact

- **Developers**: Gengrui Zhu, Yi Chen
- **Commercial licensing**: gengruizhu@outlook.com, cy150868@163.com
- **Issues**: please use the GitHub Issues page for bug reports and questions.
