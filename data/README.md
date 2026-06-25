# `data/` — input files

This directory is where GraphAllele expects its input genomes. It ships empty
on purpose (genome assemblies and annotations are large and dataset-specific).

Place the following files here, or point the run script / CLI at their absolute
paths instead.

## Required

| File | Description |
| --- | --- |
| `genome.gff3` | Whole-genome GFF3 annotation of the polyploid being analysed. |
| `genome.fasta` | Whole-genome FASTA assembly matching `genome.gff3`. |
| `reference.gff3` | Reference-genome GFF3 (used to annotate `Ref_Gene` / `Ref_Locus`). |
| `reference.cds` | Reference-genome **nucleotide** CDS FASTA (TBLASTN subject). |

## Optional

| File | Description |
| --- | --- |
| `Orthogroups.tsv` | Pre-computed OrthoFinder `Orthogroups.tsv`. Use with `-og`. If you instead pass `--auto_og`, GraphAllele runs OrthoFinder per chromosome group and you do not need this file. |

## Chromosome naming

Chromosomes **must** follow the `Chr{NUM}{SUFFIX}` convention, where `SUFFIX`
is the haplotype/subgenome letter, e.g. `Chr1A`, `Chr1B`, ..., `Chr1K`. The
`--sub_list` argument must list exactly those suffix letters.

The file names above are only examples — any names work as long as you update
the paths in `workflow/run_GraphAllele.sh` (or the `-g/-f/-ref_g/-ref_f`
arguments) accordingly.
