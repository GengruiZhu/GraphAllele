#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import argparse
import pandas as pd
import networkx as nx
from collections import defaultdict
import re

def safe_extract_id(raw_id):
    """Strip database prefixes and trailing transcript suffixes, keeping the core ID."""
    cleaned = str(raw_id).replace('gene:', '').replace('gene-', '').replace('transcript:', '').strip()
    return re.sub(r'\.[tT]?\d+$', '', cleaned)

def parse_gff(gff_file):
    """Parse GFF3 to extract gene physical coordinates for distance filtering."""
    gene_pos = {}
    chrom_order = defaultdict(list)

    if not os.path.exists(gff_file):
        print(f"[ERROR] GFF file not found: {gff_file}")
        exit(1)

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue

            info_field = cols[8]
            raw_id = None
            for field in info_field.split(';'):
                if field.startswith('ID='):
                    raw_id = field.split('=', 1)[1]
                    break

            if raw_id:
                gene_id = safe_extract_id(raw_id)
                chrom_name = cols[0]
                chrom_order[chrom_name].append(gene_id)
                gene_pos[gene_id] = (chrom_name, len(chrom_order[chrom_name]) - 1)

    print(f"[INFO] Parsed {len(gene_pos)} genes from GFF for clustering.")
    return gene_pos

def load_tandem(tandem_file):
    """Load tandem repeat blacklist."""
    blacklist = set()
    if not tandem_file or not os.path.exists(tandem_file):
        print("[WARNING] No tandem blacklist provided. Proceeding without tandem filtration.")
        return blacklist

    with open(tandem_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            genes = line.strip().split(",")
            blacklist.update([safe_extract_id(g) for g in genes])

    print(f"[INFO] Loaded {len(blacklist)} tandem duplicated genes into the blacklist.")
    return blacklist

def run_clustering(gff_file, jcvi_dir, tandem_file, output_file, min_chroms, max_dist, sub_list, chr_num):
    gene_pos = parse_gff(gff_file)
    blacklist = load_tandem(tandem_file)
    G = nx.Graph()

    anchor_files = glob.glob(os.path.join(jcvi_dir, "*.anchors"))
    if not anchor_files:
        print(f"[WARNING] No .anchors files found in {jcvi_dir}!")

    valid_edges = 0
    for fn in anchor_files:
        with open(fn, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue

                g1 = safe_extract_id(parts[0])
                g2 = safe_extract_id(parts[1])

                if g1 in gene_pos and g2 in gene_pos:
                    c1, i1 = gene_pos[g1]
                    c2, i2 = gene_pos[g2]

                    if c1 != c2:
                        G.add_edge(g1, g2)
                        valid_edges += 1
                    else:
                        if abs(i1 - i2) <= max_dist:
                            if g1 not in blacklist and g2 not in blacklist:
                                G.add_edge(g1, g2)
                                valid_edges += 1

    print(f"[INFO] Built graph with {G.number_of_nodes()} nodes and {valid_edges} edges.")

    subs = [s.strip() for s in sub_list.split(',') if s.strip()]
    standard_cols = [f"Chr{int(chr_num):02d}{s}" for s in subs]

    results = []
    for component in nx.connected_components(G):
        unique_chroms = set(gene_pos[g][0] for g in component)
        if len(unique_chroms) >= min_chroms:
            row_map = {col: [] for col in standard_cols}
            row_map["Unplaced_or_Translocated"] = []

            for gene in component:
                chrom_name = gene_pos[gene][0]
                placed = False
                expected_prefix_v1 = f"Chr{int(chr_num)}"
                expected_prefix_v2 = f"Chr{int(chr_num):02d}"

                for s, col_name in zip(subs, standard_cols):
                    if chrom_name.endswith(s) and (expected_prefix_v1 in chrom_name or expected_prefix_v2 in chrom_name):
                        row_map[col_name].append(gene)
                        placed = True
                        break

                if not placed:
                    row_map["Unplaced_or_Translocated"].append(gene)

            final_row = {}
            for k, v in row_map.items():
                final_row[k] = ",".join(v) if v else "NA"
            results.append(final_row)

    df = pd.DataFrame(results, columns=standard_cols + ["Unplaced_or_Translocated"])
    df.index = [f"Cluster_{int(chr_num):02d}_{i:05d}" for i in range(1, len(results) + 1)]
    df.index.name = "ClusterID"
    df.to_csv(output_file, sep='\t')
    print(f"[SUCCESS] Saved {len(results)} clusters to {output_file}.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Graph-based clustering of syntenic anchors into allele sets.")
    parser.add_argument('--gff', required=True, help="Merged GFF file for the chromosome group")
    parser.add_argument('--jcvi_dir', required=True, help="Directory containing .anchors files")
    parser.add_argument('--tandem', help="Path to tandem array blacklist file")
    parser.add_argument('--output', required=True, help="Output TSV file path")
    parser.add_argument('--min_chroms', type=int, default=3, help="Minimum spanning chromosomes to retain cluster")
    parser.add_argument('--max_gene_distance', type=int, default=30, help="Max index distance for intra-chromosomal edges")
    parser.add_argument('--sub_list', required=True, help="Comma-separated subgenomes (e.g., A,B,C)")
    parser.add_argument('--chr_num', required=True, help="Chromosome group number (e.g., 1)")

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    run_clustering(args.gff, args.jcvi_dir, args.tandem, args.output,
                   args.min_chroms, args.max_gene_distance, args.sub_list, args.chr_num)
