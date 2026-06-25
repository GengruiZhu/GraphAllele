#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import argparse
import logging
import sys
import networkx as nx
from collections import defaultdict
import re

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def safe_extract_id(raw_id):
    """Strip database prefixes and trailing transcript suffixes, keeping the core ID."""
    cleaned = str(raw_id).replace('gene:', '').replace('gene-', '').replace('transcript:', '').strip()
    return re.sub(r'\.[tT]?\d+$', '', cleaned)

def parse_gff_order(gff_file):
    """Index every gene by (chromosome, positional order) from a GFF3 file."""
    gene_order = {}
    chrom_genes = defaultdict(list)

    if not os.path.exists(gff_file):
        logging.error(f"GFF file not found: {gff_file}")
        sys.exit(1)

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom = parts[0]
            info = parts[8]

            gene_id = None
            for field in info.split(';'):
                if field.startswith('ID='):
                    raw_id = field.split('=', 1)[1]
                    gene_id = safe_extract_id(raw_id)
                    break

            if gene_id:
                chrom_genes[chrom].append(gene_id)

    for chrom in chrom_genes:
        for index, g_id in enumerate(chrom_genes[chrom]):
            gene_order[g_id] = (chrom, index)

    logging.info(f"Successfully indexed {len(gene_order)} genes from {len(chrom_genes)} chromosomes.")
    return gene_order

def run_self_blast(pep_file, out_blast, threads):
    """Build a clean protein DB and run an all-vs-all self-BLASTP."""
    # Ultimate aggressive FASTA cleaning: create an absolutely pristine temporary file.
    clean_pep = pep_file + ".final_clean"
    logging.info(f"Aggressively cleaning FASTA: {pep_file} -> {clean_pep}")

    with open(pep_file, 'r') as fin, open(clean_pep, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            if line.startswith('>'):
                fout.write(line + '\n')  # Retain header lines
            else:
                # Allow only letters A-Z and a-z; discard all other characters (spaces, digits, -, *, .)
                pure_seq = re.sub(r'[^A-Za-z]', '', line)
                if pure_seq:
                    fout.write(pure_seq + '\n')

    # Redirect subsequent operations to this "ultimately clean" file
    pep_file = clean_pep

    db_name = pep_file + ".db"
    logging.info(f"Building BLAST database for {pep_file}...")
    subprocess.run([
        'makeblastdb', '-in', pep_file, '-dbtype', 'prot', '-out', db_name
    ], check=True, stdout=subprocess.DEVNULL)

    logging.info(f"Running self-BLASTP with {threads} threads...")
    subprocess.run([
        'blastp', '-query', pep_file, '-db', db_name,
        '-outfmt', '6 qseqid sseqid pident evalue qcovs',
        '-evalue', '1e-10', '-num_threads', str(threads), '-out', out_blast
    ], check=True)

def identify_tandems(blast_file, gene_order, max_dist, min_ident):
    """Build a graph linking neighbouring high-identity homologs into tandem arrays."""
    G = nx.Graph()
    logging.info(f"Parsing BLAST results for tandem identification (Max Dist: {max_dist})...")

    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5: continue

            qseq = safe_extract_id(parts[0])
            sseq = safe_extract_id(parts[1])
            identity = float(parts[2])
            qcovs = float(parts[4])  # query coverage

            if qseq == sseq: continue

            if identity < min_ident or qcovs < 50.0:
                continue

            if qseq in gene_order and sseq in gene_order:
                chrom_q, pos_q = gene_order[qseq]
                chrom_s, pos_s = gene_order[sseq]

                if chrom_q == chrom_s and abs(pos_q - pos_s) <= max_dist:
                    G.add_edge(qseq, sseq)

    return G

def main():
    parser = argparse.ArgumentParser(description="Tandem Duplication Identification with Global ID Sync.")
    parser.add_argument('--gff', required=True, help='GFF3 file')
    parser.add_argument('--pep', required=True, help='Protein FASTA file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--prefix', default='species', help='Output prefix')
    parser.add_argument('--threads', type=int, default=8, help='Threads for BLAST')
    parser.add_argument('--max_dist', type=int, default=5, help='Max gene index distance (Default: 5)')
    parser.add_argument('--min_ident', type=float, default=50.0, help='Min identity (Default: 50.0)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    blast_out = os.path.join(args.outdir, f"{args.prefix}.blast")
    tandem_out = os.path.join(args.outdir, f"{args.prefix}.tandem")

    if os.path.exists(tandem_out) and os.path.getsize(tandem_out) > 0:
        logging.info("Tandem file exists, skipping.")
        return

    gene_order = parse_gff_order(args.gff)
    run_self_blast(args.pep, blast_out, args.threads)
    tandem_graph = identify_tandems(blast_out, gene_order, args.max_dist, args.min_ident)

    with open(tandem_out, 'w') as fout:
        for component in nx.connected_components(tandem_graph):
            fout.write(",".join(sorted(list(component))) + "\n")

    logging.info(f"Found {tandem_graph.number_of_nodes()} tandem duplicated genes. Results saved to {tandem_out}")

if __name__ == '__main__':
    main()
