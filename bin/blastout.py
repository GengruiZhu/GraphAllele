#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import argparse
import tempfile
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import re

def safe_extract_id(raw_id):
    """
    Safely strips database prefixes (gene:, transcript:) and
    trailing transcript suffixes (e.g., .1, .2, .t1, .T01),
    preserving the core biological ID and internal dots.
    """
    cleaned = str(raw_id).replace('gene:', '').replace('gene-', '').replace('transcript:', '').strip()
    return re.sub(r'\.[tT]?\d+$', '', cleaned)

def run_blast(query, subject, out, threads, temp_dir):
    """
    Execute BLASTP for protein sequence alignment within a secure sandbox.
    The subject FASTA is aggressively cleaned and written into the temp
    directory to prevent polluting the user's original reference data
    directory with makeblastdb artifacts.
    """
    local_subject = os.path.join(temp_dir, "subject.fa")

    # Read the original subject FASTA and write a perfectly clean version to the temp dir
    with open(subject, 'r') as src, open(local_subject, 'w') as dst:
        for line in src:
            line_strip = line.strip()
            if not line_strip:
                continue  # Skip empty lines
            if line_strip.startswith('>'):
                dst.write(line_strip + '\n')  # Preserve headers exactly as they are
            else:
                # Aggressively remove any non-alphabet character (*, -, space, dots, etc.)
                clean_seq = re.sub(r'[^A-Za-z]', '', line_strip)
                if clean_seq:
                    dst.write(clean_seq + '\n')

    subprocess.run(['makeblastdb', '-in', local_subject, '-dbtype', 'prot'],
                   check=True, stdout=subprocess.DEVNULL)

    subprocess.run(['blastp', '-query', query, '-db', local_subject,
                    '-outfmt', '6 qseqid sseqid pident qcovs',
                    '-num_threads', str(threads), '-out', out, '-evalue', '1e-10'],
                   check=True)

def main():
    parser = argparse.ArgumentParser(description="True Sequence Homology Expansion via BLASTP.")
    parser.add_argument('--allele_file', required=True, help='Input verified clustered TSV file.')
    parser.add_argument('--fasta', required=True, help='Merged Protein FASTA file.')
    parser.add_argument('--out_prefix', required=True, help='Prefix for output files.')
    parser.add_argument('--threads', type=int, default=10, help='Number of threads for BLAST.')
    parser.add_argument('--identity', type=float, default=90.0,
                        help='Minimum sequence identity to salvage translocated genes (Default: 90.0).')
    args = parser.parse_args()

    df = pd.read_csv(args.allele_file, sep='\t')
    # Resolve FutureWarning type conflicts: force the target column to string type
    # and replace original NaN values (recognized as float64) with empty strings.
    target_col = 'Unplaced_or_Translocated'
    if target_col in df.columns:
        df[target_col] = df[target_col].astype(str).replace(['nan', 'NaN', 'NA'], '')
    else:
        df[target_col] = ""

    # Safely identify structural columns
    allele_cols = [c for c in df.columns if c.startswith('Chr') or c == 'Unplaced_or_Translocated']

    # Map clustered genes to their row index for later expansion insertion
    gene_to_row_idx = {}
    all_clustered_genes = set()

    for idx, row in df.iterrows():
        for c in allele_cols:
            if c not in df.columns or pd.isna(row[c]) or str(row[c]).strip() in ('NA', 'nan', ''):
                continue
            for g in str(row[c]).split(','):
                clean_id = safe_extract_id(g)
                if clean_id:
                    gene_to_row_idx[clean_id] = idx
                    all_clustered_genes.add(clean_id)

    # Use a secure temporary directory that auto-cleans everything (query, db artifacts, blast output)
    with tempfile.TemporaryDirectory() as tmp:
        query_fa = os.path.join(tmp, "query.fa")
        blast_res = os.path.join(tmp, "blast_results.tmp")

        extracted_count = 0
        with open(query_fa, "w", encoding="utf-8") as out:
            for rec in SeqIO.parse(args.fasta, "fasta"):
                seq_id = safe_extract_id(rec.id)
                if seq_id in all_clustered_genes:
                    # Write sequence with the cleaned ID for downstream perfect matching
                    rec.id = seq_id
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
                    extracted_count += 1

        print(f"[INFO] Extracted {extracted_count} target sequences for homology expansion.")

        if extracted_count == 0:
            print("[WARNING] Zero sequences extracted. Skipping BLAST. Please check FASTA headers.")
            df.to_csv(args.out_prefix + "_expanded.tsv", sep='\t', index=False)
            return

        print("[INFO] Building database and running BLASTP refinement...")
        run_blast(query_fa, args.fasta, blast_res, args.threads, tmp)

        # Parse BLAST results to identify high-identity unclustered paralogs/alleles
        salvaged_genes_count = 0
        expansion_map = defaultdict(set)

        if os.path.exists(blast_res):
            best_salvage = {}
            with open(blast_res, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < 4: continue

                    qseq, sseq = safe_extract_id(parts[0]), safe_extract_id(parts[1])
                    pident = float(parts[2])
                    qcovs = float(parts[3])
                    # Ignore self-hits
                    if qseq == sseq: continue

                    # The core expansion logic:
                    # If an unclustered gene hits a clustered gene with high identity...
                    if pident >= args.identity and qcovs >= 50.0 and sseq not in all_clustered_genes:
                        if qseq in gene_to_row_idx:
                            target_row = gene_to_row_idx[qseq]
                            if sseq not in best_salvage or pident > best_salvage[sseq]['pident']:
                                best_salvage[sseq] = {'pident': pident, 'target_row': target_row}

            for sseq, info in best_salvage.items():
                expansion_map[info['target_row']].add(sseq)

        # Update the DataFrame with salvaged genes
        if expansion_map:
            for idx, new_genes in expansion_map.items():
                existing_unplaced = str(df.at[idx, 'Unplaced_or_Translocated'])
                salvaged_genes_count += len(new_genes)

                if existing_unplaced in ('NA', 'nan', ''):
                    df.at[idx, 'Unplaced_or_Translocated'] = ",".join(new_genes)
                else:
                    # Merge existing unplaced genes with newly salvaged ones, removing duplicates
                    combined = set(safe_extract_id(g) for g in existing_unplaced.split(',')) | new_genes
                    df.at[idx, 'Unplaced_or_Translocated'] = ",".join(sorted(list(combined)))

            print(f"[SUCCESS] Salvaged {salvaged_genes_count} missing homologous genes into the matrix!")
        else:
            print("[INFO] No novel translocated/unplaced genes passed the identity threshold for expansion.")

        df.to_csv(args.out_prefix + "_expanded.tsv", sep='\t', index=False)
        print(f"[INFO] Expansion finished. Output generated: {args.out_prefix}_expanded.tsv")

if __name__ == '__main__':
    main()
