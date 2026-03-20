#!/bin/bash
#PBS -N GraphAllele
#PBS -l nodes=1:ppn=20
#PBS -q bigfat
#PBS -l mem=
#PBS -o Allele.log
#PBS -j oe

cd $PBS_O_WORKDIR

date -R

source ~/miniconda3/bin/activate

source activate polyalleler

python GraphAllele.py \
  -g ../data/ROC22.gff3 \
  -f ../data/ROC22.fasta \
  -ref_g ../data/Eru.gff3 \
  -ref_f ../data/Eru.cds \
  -og ./Orthogroups.tsv \
  -s 1 -e 10 \
  -t 20 \
  --sub_list A,B,C,D,E,F,G,H,I,J,K,L,M,N
