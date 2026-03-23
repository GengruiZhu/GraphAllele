#!/bin/bash
#PBS -N GraphAllele_v2
#PBS -l nodes=1:ppn=30
#PBS -q bigfat
#PBS -l mem=
#PBS -o Allele.log
#PBS -j oe

cd $PBS_O_WORKDIR

date -R

source ~/miniconda3/bin/activate

source activate polyalleler

python GraphAllele_v2.py \
  -g ../data/ROC22.gff3 \
  -f ../data/ROC22.fasta \
  -ref_g ../data/Eru.gff3 \
  -ref_f ../data/Eru.cds \
  --auto_og \
  -s 1 -e 10 \
  -t 30 \
  --sub_list A,B,C,D,E,F,G,H,I,J,K,L,M,N
