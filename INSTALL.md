\# Installation Guide



This document provides detailed installation instructions for GraphAllele and all its dependencies.



\---



\## System Requirements



| Component | Minimum | Recommended |

|---|---|---|

| OS | Linux (any major distro) | CentOS 7 / Ubuntu 20.04+ |

| CPU | 4 cores | 20+ cores |

| RAM | 32 GB | 128 GB+ (for large polyploid genomes) |

| Disk | 100 GB free | 500 GB+ |

| Python | 3.8 | 3.9 or 3.10 |

| Conda | Any | Miniconda3 or Anaconda3 |



> \*\*Note\*\*: For highly polyploid genomes (e.g., sugarcane with 100+ chromosomes), RAM usage can exceed 64 GB during tandem identification and BLAST steps. High-memory compute nodes are strongly recommended.



\---



\## Step 1: Install Conda (if not already installed)



```bash

\# Download Miniconda3 installer

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86\_64.sh



\# Run the installer

bash Miniconda3-latest-Linux-x86\_64.sh



\# Reload shell

source \~/.bashrc

```



\---



\## Step 2: Clone GraphAllele



```bash

git clone https://github.com/GengruiZhu/GraphAllele.git

cd GraphAllele

```



\---



\## Step 3: Create the Conda Environment



\### Option A: Use the provided environment file (Recommended)



```bash

conda env create -f environment.yml

conda activate polyalleler

```



\### Option B: Manual step-by-step creation



```bash

\# Create environment

conda create -n polyalleler python=3.9 -y

conda activate polyalleler



\# Install Python packages

pip install biopython pandas networkx



\# Install bioinformatics tools via bioconda

conda install -c bioconda -c conda-forge blast gffread last jcvi -y



\# Optional: OrthoFinder (only needed if using --auto\_og)

conda install -c bioconda orthofinder -y

```



\---



\## Step 4: Verify Installation



Run the following to confirm all dependencies are available:



```bash

\# Activate environment first

conda activate polyalleler



\# Check Python packages

python -c "from Bio import SeqIO; import pandas; import networkx; print('Python packages: OK')"



\# Check external tools

blastp -version

makeblastdb -version

tblastn -version

gffread --version

python -m jcvi.compara.catalog --help 2>\&1 | head -3

```



If all commands print version info without errors, installation is complete.



\---



\## Step 5: Test with Example Data (Optional)



```bash

cd GraphAllele/workflow



\# Test run with provided example data (Chr5 only)

python GraphAllele.py \\

&#x20; -g ../data/ROC22.gff3 \\

&#x20; -f ../data/ROC22.fasta \\

&#x20; -ref\_g ../data/Eru.gff3 \\

&#x20; -ref\_f ../data/Eru.cds \\

&#x20; -og ./Orthogroups.tsv \\

&#x20; -s 5 -e 5 \\

&#x20; -t 4 \\

&#x20; --sub\_list A,B,C,D,E

```



\---



\## Troubleshooting



\### `ModuleNotFoundError: No module named 'Bio'`



```bash

conda activate polyalleler

pip install biopython

```



\### `blastp: command not found`



```bash

conda activate polyalleler

conda install -c bioconda blast -y

```



\### `gffread: command not found`



```bash

conda activate polyalleler

conda install -c bioconda gffread -y

```



\### JCVI/MCScan errors during Step 3



Ensure you are running from within the correct working directory. GraphAllele handles this automatically via `cwd=` in subprocess calls. If issues persist, check that `.cds` and `.bed` files are correctly linked in the `03.jcvi/` subdirectory.



\### OrthoFinder not found (when using `--auto\_og`)



```bash

conda activate polyalleler

conda install -c bioconda orthofinder -y

```



\---



\## PBS Job Submission



Edit `workflow/run\_GraphAllele.sh` to set the correct `mem=` and `ppn=` values for your cluster:



```bash

\#PBS -l nodes=1:ppn=20

\#PBS -l mem=128gb

```



\---



\## HPC / Cluster Notes



\- The pipeline processes chromosome groups \*\*sequentially\*\* (one by one) to protect HPC NFS I/O and ensure stable logging. Set `-t` to the number of threads available per node.

\- For PBS clusters, submit via `qsub workflow/run\_GraphAllele.sh`.

\- For SLURM clusters, adapt the job header accordingly:



```bash

\#!/bin/bash

\#SBATCH -J GraphAllele

\#SBATCH -n 20

\#SBATCH --mem=128G

\#SBATCH -o Allele.log



source activate polyalleler

cd /path/to/GraphAllele/workflow



python GraphAllele.py \\

&#x20; -g ../data/genome.gff3 \\

&#x20; -f ../data/genome.fasta \\

&#x20; -ref\_g ../data/reference.gff3 \\

&#x20; -ref\_f ../data/reference.cds \\

&#x20; --auto\_og \\

&#x20; -s 1 -e 10 \\

&#x20; -t 20 \\

&#x20; --sub\_list A,B,C,D,E,F,G,H,I,J,K,L,M,N

```

