#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import argparse
import tempfile
import pandas as pd
import re
from collections import Counter

def safe_extract_id(raw_id):
    """
    Safely strips database prefixes (gene:, transcript:) and
    trailing transcript suffixes (e.g., .1, .2, .t1, .T01),
    preserving the core biological ID and internal dots.
    """
    cleaned = str(raw_id).replace('gene:', '').replace('gene-', '').replace('transcript:', '').strip()
    return re.sub(r'\.[tT]?\d+$', '', cleaned)

def parse_gff_locus(gff_file):
    """Parse reference GFF3 to extract accurate physical loci."""
    locus_map = {}

    if not os.path.exists(gff_file):
        print(f"[ERROR] Reference GFF file not found: {gff_file}")
        exit(1)

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) > 8 and parts[2] == "gene":
                try:
                    # Extract ID attribute
                    info = parts[8]
                    gene_id = None
                    for field in info.split(';'):
                        if field.startswith('ID='):
                            gene_id = field.split('=', 1)[1]
                            break

                    if gene_id:
                        clean_id = safe_extract_id(gene_id)
                        locus = f"{parts[0]}:{parts[3]}-{parts[4]}({parts[6]})"
                        locus_map[clean_id] = locus
                except IndexError:
                    continue

    print(f"[INFO] Parsed {len(locus_map)} distinct loci from reference GFF.")
    return locus_map

def run_tblastn(query_fasta, subject_fasta, blast_out, threads, temp_dir):
    """
    Executes TBLASTN securely within a sandbox to prevent database
    artifacts (.nin, .nhr) from polluting the shared reference directory.
    """
    local_subject = os.path.join(temp_dir, "ref_subject.fa")

    # Symlink or copy the reference into the sandbox
    try:
        os.symlink(os.path.abspath(subject_fasta), local_subject)
    except OSError:
        import shutil
        shutil.copy2(subject_fasta, local_subject)

    subprocess.run(['makeblastdb', '-in', local_subject, '-dbtype', 'nucl'],
                   check=True, stdout=subprocess.DEVNULL)

    # Run TBLASTN with an explicit E-value to prevent spurious short alignments
    subprocess.run(['tblastn', '-query', query_fasta, '-db', local_subject,
                    '-outfmt', '6 qseqid sseqid pident',
                    '-num_threads', str(threads), '-out', blast_out, '-evalue', '1e-10'],
                   check=True, stderr=subprocess.PIPE)

def main():
    parser = argparse.ArgumentParser(description="Map haplotype clusters back to reference genome locus.")
    parser.add_argument('--allele_file', required=True, help='Input expanded cluster TSV file.')
    parser.add_argument('--ref_gff', required=True, help='Reference GFF3 file.')
    parser.add_argument('--ref_cds', required=True, help='Reference Nucleotide CDS FASTA.')
    parser.add_argument('--hap_cds', required=True, help='Haplotype Protein CDS FASTA (m_pep).')
    parser.add_argument('--output', required=True, help='Output standardized TSV file.')
    parser.add_argument('--sub_list', required=True, help='Comma-separated subgenome list (e.g., A,B,C).')
    parser.add_argument('--chr_num', required=True, help='Base chromosome number.')
    parser.add_argument('--identity', type=float, default=80.0, help='Minimum identity for TBLASTN mapping.')
    parser.add_argument('--threads', type=int, default=10, help='Number of threads for TBLASTN.')
    args = parser.parse_args()

    print("[INFO] Parsing reference GFF loci...")
    locus_map = parse_gff_locus(args.ref_gff)

    # Use a secure sandbox for the BLAST mapping phase
    with tempfile.TemporaryDirectory() as tmp:
        blast_tmp = os.path.join(tmp, "tblastn_results.tmp")

        print(f"[INFO] Running TBLASTN to map haplotype proteins to reference nucleotides...")
        run_tblastn(args.hap_cds, args.ref_cds, blast_tmp, args.threads, tmp)

        # Dictionary to hold the strongest biological hit for a query
        # Structure: best_hits[query_id] = (subject_id, pident)
        best_hits = {}

        if os.path.exists(blast_tmp) and os.path.getsize(blast_tmp) > 0:
            b_df = pd.read_csv(blast_tmp, sep='\t', header=None)

            # Filter by minimum identity threshold
            b_df = b_df[b_df[2] >= args.identity]

            for _, r in b_df.iterrows():
                q_id = safe_extract_id(r[0])
                s_id = safe_extract_id(r[1])
                pident = float(r[2])

                # Maintain the highest identity mapping for each query gene
                if q_id not in best_hits or pident > best_hits[q_id][1]:
                    best_hits[q_id] = (s_id, pident)
        else:
            print("[WARNING] TBLASTN produced no valid alignments. The cluster matrix will lack reference coordinates.")

    df = pd.read_csv(args.allele_file, sep='\t')

    # Identify all structural columns to prevent data loss
    allele_cols = [c for c in df.columns if c.startswith('Chr') or c == 'Unplaced_or_Translocated']
    final_rows = []

    print("[INFO] Projecting best reference coordinates onto haplotype clusters...")

    for index, row in df.iterrows():
        cluster_id = row['ClusterID'] if 'ClusterID' in row else f"Cluster_{len(final_rows):05d}"

        cluster_ref_votes = []

        for col in allele_cols:
            if col in row and pd.notna(row[col]) and str(row[col]).strip() not in ('NA', 'nan', ''):
                cell_content = str(row[col])

                for g in cell_content.split(','):
                    gene = safe_extract_id(g)

                    if gene in best_hits:
                        candidate_ref_id, _ = best_hits[gene]
                        cluster_ref_votes.append(candidate_ref_id)

        if cluster_ref_votes:
            most_common_ref_id = Counter(cluster_ref_votes).most_common(1)[0][0]
            best_ref_id = most_common_ref_id
            best_ref_locus = locus_map.get(most_common_ref_id, "NA")
        else:
            best_ref_id = "NA"
            best_ref_locus = "NA"

        new_row = [cluster_id, best_ref_id, best_ref_locus] + [str(row.get(col, "NA")) for col in allele_cols]
        final_rows.append(new_row)

    out_cols = ['ClusterID', 'Ref_Gene', 'Ref_Locus'] + allele_cols
    final_df = pd.DataFrame(final_rows, columns=out_cols)
    final_df.to_csv(args.output, sep='\t', index=False)

    print(f"[SUCCESS] Final standardized mapping table generated: {os.path.abspath(args.output)}")

if __name__ == '__main__':
    main()
