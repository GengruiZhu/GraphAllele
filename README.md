\# GraphAllele



A graph-constrained pipeline for constructing standardized allele matrices in polyploid genomes.



\---



\## Table of Contents



\- \[Overview](#overview)

\- \[Features](#features)

\- \[Pipeline Architecture](#pipeline-architecture)

\- \[Installation](#installation)

\- \[Quick Start](#quick-start)

\- \[Usage](#usage)

\- \[Output Description](#output-description)

\- \[Dependencies](#dependencies)

\- \[Changelog](#changelog)

\- \[License](#license)

\- \[Contact](#contact)



\---



\## Overview



GraphAllele integrates synteny-based graph clustering, tandem duplication filtering, orthogroup verification, and reference genome calibration into a unified, resumable pipeline for allele identification in complex polyploid genomes.



\*\*Chromosome naming convention required\*\*: `Chr{NUM}{SUFFIX}`, e.g., `Chr1A`, `Chr1B`, ..., `Chr1K`.



\---



\## Features



\- \*\*Breakpoint-resumable\*\*: each step is checkpointed and skipped automatically on re-run; if the pipeline fails at any step, simply re-run the same command and it will resume from the last incomplete step.

\- \*\*Sequential processing\*\*: chromosome groups are processed one by one to protect HPC NFS I/O and ensure stable, readable logging.

\- \*\*Graph-constrained clustering\*\*: NetworkX-based extraction of syntenic connected components with configurable distance thresholds.

\- \*\*Tandem duplication filtering\*\*: self-BLASTP + graph-based tandem array blacklisting with adjustable gene distance.

\- \*\*Intra-group OrthoFinder with early termination\*\*: OrthoFinder is run independently per chromosome group; the process is killed immediately after `Orthogroups.tsv` is written to disk, skipping the MSA and tree-building phases to reduce runtime.

\- \*\*Reference calibration\*\*: TBLASTN-based anchoring to an external reference genome for cross-species annotation.

\- \*\*Standardized output\*\*: fixed-column allele matrix (A–N subgenomes) with globally unique cluster IDs.

\- \*\*Automated JCVI ortholog analysis\*\*: built-in shell script for pairwise synteny comparisons with auto-discovery of species from input files.

\- \*\*Global matrix normalization\*\*: post-processing module that produces both a cleaned global matrix and a k-mer–formatted cluster file for downstream analysis.



\---



\## Pipeline Architecture



```

Input: Whole-genome GFF3 + FASTA + Reference CDS/GFF + Orthogroups.tsv (or --auto\_og)

&#x20;        |

&#x20;        v

&#x20; Step 1: Prepare Data (01.prepare/)

&#x20;         Split GFF/FASTA by chromosome group

&#x20;         Extract CDS, PEP, BED per haplotype (gffread)

&#x20;        |

&#x20;        v

&#x20; Step 1.5 (optional): Intra-Group OrthoFinder (--auto\_og)

&#x20;         Run OrthoFinder per chromosome group

&#x20;         Early termination after Orthogroups.tsv is written

&#x20;        |

&#x20;        v

&#x20; Step 2: Tandem Duplication Identification (02.tandem/)

&#x20;         Self-BLASTP on merged PEP

&#x20;         Graph-based tandem array detection (--tandem\_dist)

&#x20;        |

&#x20;        v

&#x20; Step 3: JCVI Synteny Analysis (03.jcvi/)

&#x20;         Pairwise LAST alignment via run\_ortholog.sh

&#x20;         MCScan anchor generation (.anchors files)

&#x20;        |

&#x20;        v

&#x20; Step 4: Graph Clustering (04.cluster.tsv)

&#x20;         Build synteny graph from .anchors

&#x20;         Filter tandem blacklist

&#x20;         Extract connected components as candidate allele clusters

&#x20;        |

&#x20;        v

&#x20; Step 5: OrthoGroup Verification (05.verified.tsv)

&#x20;         Validate clusters against OrthoFinder orthogroups

&#x20;         Output verified and rejected tables separately

&#x20;        |

&#x20;        v

&#x20; Step 6: Sequence Homology Expansion (06.expanded/)

&#x20;         BLASTP-based refinement of allele assignments

&#x20;        |

&#x20;        v

&#x20; Step 7: Reference Calibration (07.FINAL\_ALLELE.tsv)

&#x20;         TBLASTN anchoring to reference CDS

&#x20;         Annotate each cluster with Ref\_Gene and Ref\_Locus

&#x20;        |

&#x20;        v

&#x20; Post-processing: Global Merge \& Normalization

&#x20;         Merge all chromosome groups into PolyAlleler\_Global\_Matrix.tsv

&#x20;         Normalize into PolyAlleler\_Global\_Matrix\_Cleaned.tsv (fixed A–N columns)

&#x20;         Generate my\_clusters.tsv (k-mer format)

```



\---



\## Installation



\### Conda (Recommended)



```bash

git clone https://github.com/GengruiZhu/GraphAllele.git

cd GraphAllele



conda env create -f environment.yml

conda activate polyalleler

```



\### Manual Installation



See \[INSTALL.md](INSTALL.md) for step-by-step instructions.



\---



\## Quick Start



\### Basic run with a pre-computed Orthogroups file



```bash

cd GraphAllele/workflow



python GraphAllele.py \\

&#x20; -g /path/to/genome.gff3 \\

&#x20; -f /path/to/genome.fasta \\

&#x20; -ref\_g /path/to/reference.gff3 \\

&#x20; -ref\_f /path/to/reference.cds \\

&#x20; -og /path/to/Orthogroups.tsv \\

&#x20; -s 1 -e 5 \\

&#x20; -t 20 \\

&#x20; --sub\_list A,B,C,D,E

```



\### Run with automatic OrthoFinder (no pre-computed file needed)



```bash

python GraphAllele.py \\

&#x20; -g /path/to/genome.gff3 \\

&#x20; -f /path/to/genome.fasta \\

&#x20; -ref\_g /path/to/reference.gff3 \\

&#x20; -ref\_f /path/to/reference.cds \\

&#x20; --auto\_og \\

&#x20; -s 1 -e 10 \\

&#x20; -t 30 \\

&#x20; --min\_c 2 \\

&#x20; --cluster\_dist 100 \\

&#x20; --sub\_list A,B,C,D,E,F,G,H,I,J,K,L,M,N

```



\### PBS Job Submission



```bash

\# Edit workflow/run\_GraphAllele.sh to set your paths and resource requirements, then:

qsub workflow/run\_GraphAllele.sh

```



\*\*Example PBS script\*\* (`workflow/run\_GraphAllele.sh`):



```bash

\#!/bin/bash

\#PBS -N GraphAllele\_ZG

\#PBS -l nodes=1:ppn=30

\#PBS -q comput

\#PBS -l mem=

\#PBS -o Allele.log

\#PBS -j oe



cd $PBS\_O\_WORKDIR

date -R

source \~/miniconda3/bin/activate

source activate polyalleler



python GraphAllele.py \\

&#x20; -g ../data/ZG.gff3 \\

&#x20; -f ../data/ZG.fasta \\

&#x20; -ref\_g ../data/Eru.gff3 \\

&#x20; -ref\_f ../data/Eru.cds \\

&#x20; --auto\_og \\

&#x20; -s 1 -e 10 \\

&#x20; -t 30 \\

&#x20; -o ZG\_100 \\

&#x20; --min\_c 2 \\

&#x20; --cluster\_dist 100 \\

&#x20; --sub\_list A,B,C,D,E,F,G,H,I,J,K,L,M,N

```



\---



\## Usage



```

python GraphAllele.py \[OPTIONS]

```



\### Required Arguments



| Argument | Description |

|---|---|

| `-g / --gff` | Whole-genome GFF3 annotation file. Must use `Chr{NUM}{SUFFIX}` naming convention. |

| `-f / --fasta` | Whole-genome FASTA sequence file corresponding to the GFF3. |

| `-ref\_g / --ref\_gff` | Reference genome GFF3 annotation (used for calibration in Step 7). |

| `-ref\_f / --ref\_cds` | Reference genome CDS FASTA (used for TBLASTN in Step 7). |



\### Optional Arguments



| Argument | Default | Description |

|---|---|---|

| `-og / --orthogroups` | `None` | Path to a pre-computed OrthoFinder `Orthogroups.tsv`. If provided, the pipeline uses this file for orthogroup verification (Step 5). Mutually exclusive in spirit with `--auto\_og`, though both can be set (pre-computed file takes priority if the auto-generated one does not exist). |

| `--auto\_og` | `False` | Run OrthoFinder automatically per chromosome group using only the haplotype PEP files for that group. The OrthoFinder process is terminated immediately after `Orthogroups.tsv` is written, skipping MSA and tree-building. Recommended when no pre-computed orthogroups are available. |

| `-s / --start` | `1` | Start chromosome group number. For example, `-s 1` means beginning from `Group\_Chr01`. |

| `-e / --end` | `10` | End chromosome group number. For example, `-e 10` means processing through `Group\_Chr10`. |

| `-t / --threads` | `10` | Total CPU threads available. Used by BLAST, OrthoFinder, and other parallelizable steps. |

| `-o / --outdir` | `standardized\_results` | Output directory. All intermediate and final results are stored here. |

| `--sub\_list` | `A,B,...,N` | Comma-separated haplotype suffix list. Defines the column order in the output allele matrix. Must match the suffix letters in your chromosome names. |

| `--min\_c` | `3` | Minimum number of distinct haplotype chromosomes a cluster must span to be retained. Lower values (e.g., 2) are more permissive and suitable for genomes with high gene loss. |

| `--tandem\_dist` | `5` | Maximum gene index distance for tandem duplicate detection (MCScanX convention). Two genes within this distance on the same chromosome with high sequence similarity are flagged as tandem duplicates. |

| `--cluster\_dist` | `30` | Maximum gene index distance for synteny graph clustering. Increasing this value (e.g., to 100) allows more distant syntenic gene pairs to be grouped together, which can be useful for genomes with many structural rearrangements. |



\### Parameter Guidance



\- \*\*`--sub\_list`\*\*: must exactly match the suffix letters used in your chromosome names. For example, if your chromosomes are named `Chr1A` through `Chr1K`, use `--sub\_list A,B,C,D,E,F,G,H,I,J,K`. The order of suffixes determines column order in the output matrix.

\- \*\*`--min\_c`\*\*: a cluster must span at least this many distinct haplotype chromosomes to be retained. For highly fragmented or recently diverged polyploids, consider lowering to `2`.

\- \*\*`--auto\_og`\*\*: when enabled, OrthoFinder is run once per chromosome group using only the haplotype PEP files for that group. The process is terminated automatically after `Orthogroups.tsv` is written, before MSA and tree steps. It is still recommended to supply a pre-computed `-og` file when one is available, as it tends to be more comprehensive.

\- \*\*`--cluster\_dist`\*\*: the default value of `30` works well for most plant genomes. For genomes with extensive rearrangements (e.g., sugarcane), values of `50–100` may produce better results.

\- \*\*`--tandem\_dist`\*\*: the default of `5` follows the MCScanX convention. Increase to `10` if you suspect many dispersed tandem arrays are being missed.



\---



\## Output Description



```

standardized\_results/

├── Group\_Chr01/

│   ├── 01.prepare/                    # Per-haplotype GFF, FASTA, CDS, BED, PEP

│   ├── 01.5.OrthoFinder\_Intra/        # Intra-group OrthoFinder results (--auto\_og only)

│   ├── 02.tandem/                     # Merged GFF/PEP + tandem blacklist

│   ├── 03.jcvi/                       # Pairwise .anchors files

│   ├── 04.cluster.tsv                 # Raw synteny clusters

│   ├── 05.verified.tsv                # OrthoGroup-validated clusters

│   ├── 05.verified\_rejected.tsv       # Clusters failing OG verification

│   ├── 06.expanded\_expanded.tsv       # BLAST-refined allele table

│   └── 07.FINAL\_ALLELE.tsv            # Final allele matrix with reference anchors

├── Group\_Chr02/

│   └── ...

├── PolyAlleler\_Global\_Matrix.tsv       # Genome-wide merged allele matrix (raw)

├── PolyAlleler\_Global\_Matrix\_Cleaned.tsv  # Normalized matrix with fixed A–N columns

└── my\_clusters.tsv                     # K-mer format cluster file

```



\### `07.FINAL\_ALLELE.tsv` column format



| ClusterID | Ref\_Gene | Ref\_Locus | Chr01A | Chr01B | ... |

|---|---|---|---|---|---|

| Global\_Cluster\_000001 | Gene001 | Chr1:1000-5000(+) | gene\_A | gene\_B | ... |



Each row is one allele group. `-` or `NA` indicates no gene assigned for that haplotype in this cluster.



\### `PolyAlleler\_Global\_Matrix\_Cleaned.tsv`



A post-processed version of the global matrix with standardized columns (`Allele\_A` through `Allele\_N`) and globally re-numbered cluster IDs, produced by the normalization module.



\### `my\_clusters.tsv`



A two-column file: `ClusterID` and a comma-separated list of all gene members in that cluster. Useful for downstream k-mer or expression analyses.



\---



\## Dependencies



\### Python Packages



| Package | Version |

|---|---|

| Python | >= 3.8 |

| Biopython | >= 1.79 |

| pandas | >= 1.3 |

| networkx | >= 2.6 |



\### External Tools



| Tool | Purpose |

|---|---|

| \[BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | Sequence similarity search (BLASTP, TBLASTN, makeblastdb) |

| \[gffread](https://github.com/gpertea/gffread) | CDS/protein extraction from GFF |

| \[JCVI/MCScan](https://github.com/tanghaibao/jcvi) | Synteny analysis |

| \[LAST](https://gitlab.com/mcfrith/last) | Pairwise sequence alignment |

| \[OrthoFinder](https://github.com/davidemms/OrthoFinder) \*(optional)\* | Orthogroup inference (required if using `--auto\_og`) |



\---



\## Changelog



See \[CHANGELOG.md](CHANGELOG.md) for a detailed version history.



\### v1.1 (2025-04-02)



\- Added automated JCVI ortholog shell script (`bin/run\_ortholog.sh`) with species auto-discovery

\- Added global matrix normalization module (`generate\_clean\_clusters\_auto`) producing cleaned matrix and k-mer cluster file

\- Improved global merge logic using pandas with re-numbered ClusterIDs

\- Added `--protein` flag to `prepare\_jcvi.py` for protein sequence extraction

\- Enhanced error handling with `\[FATAL ERROR]` / `\[ACTION]` messages and resume guidance

\- Updated PBS example script with real-world parameters



\### v1.0



\- Initial release



\---



\## Citation



> Manuscript in preparation. Citation will be added upon publication.



\---



\## License



This project is licensed under a \*\*Non-Commercial Research License\*\*.

Free to use for academic and research purposes only. Commercial use is strictly prohibited.

See the \[LICENSE](LICENSE) file for details.



\---



\## Contact



\- \*\*Developers\*\*: Gengrui Zhu, Yi Chen

\- \*\*Issues\*\*: please use the \[GitHub Issues](https://github.com/GengruiZhu/GraphAllele/issues) page for bug reports and questions.

