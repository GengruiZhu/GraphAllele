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

python GraphAllele_final.py \
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
