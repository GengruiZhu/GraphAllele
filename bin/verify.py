#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import os
import sys
import re
from collections import defaultdict

def safe_extract_id(raw_id):
    """
    Safely strips database prefixes (gene:, transcript:) and
    trailing transcript suffixes (e.g., .1, .2, .t1, .T01),
    preserving the core biological ID and internal dots.
    """
    cleaned = str(raw_id).replace('gene:', '').replace('gene-', '').replace('transcript:', '').strip()
    return re.sub(r'\.[tT]?\d+$', '', cleaned)

def extract_subgenome_from_id(gene_id):
    """Infer the subgenome suffix (e.g. 'A') embedded in a gene ID like Chr1A_g0001."""
    match = re.search(r'Chr\d+([a-zA-Z]+)', gene_id, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return None

def load_og(og_file):
    """Load Gene to Orthogroup mapping dictionary robustly"""
    og_map = {}
    og_to_genes = defaultdict(set)
    if not os.path.exists(og_file):
        print(f"[WARNING] Orthogroups file {og_file} not found.")
        return og_map, og_to_genes

    try:
        df = pd.read_csv(og_file, sep='\t', low_memory=False)
        for _, row in df.iterrows():
            og_id = row['Orthogroup']
            for genes in row[1:].dropna():
                if pd.isna(genes) or str(genes).strip() == "":
                    continue
                for g in str(genes).split(','):
                    clean_id = safe_extract_id(g)
                    if clean_id:
                        og_map[clean_id] = og_id
                        og_to_genes[og_id].add(clean_id)

    except Exception as e:
        print(f"[ERROR] Failed to parse Orthogroups file: {e}")
        sys.exit(1)

    print(f"[INFO] Loaded {len(og_map)} distinct gene-to-OG mappings from OrthoFinder.")
    return og_map, og_to_genes

def extract_genes_from_row(row, target_cols):
    """Safely extract all individual genes, handling comma-separated paralogs."""
    genes = []
    for c in target_cols:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip() not in ('NA', '-', '', 'nan'):
            cell_content = str(row[c])
            for g in cell_content.split(','):
                clean_gene = safe_extract_id(g)
                if clean_gene:
                    genes.append(clean_gene)
    return genes

def main():
    parser = argparse.ArgumentParser(description="Orthogonal Verification of Syntenic Clusters (With Native Rescue)")
    parser.add_argument('--cluster_file', required=True, help="Input synteny clusters TSV")
    parser.add_argument('--orthogroups', required=True, help="Input OrthoFinder TSV")
    parser.add_argument('--output_prefix', required=True, help="Prefix for verified/rejected TSVs")
    parser.add_argument('--min_support_ratio', type=float, default=0.35,
                        help='Minimum ratio of genes sharing the dominant OG. Default: 0.35')
    args = parser.parse_args()

    og_map, og_to_genes = load_og(args.orthogroups)
    if not og_map:
        print("[ERROR] Orthogroup map is empty. Aborting verification to prevent false rejections.")
        sys.exit(1)

    try:
        df = pd.read_csv(args.cluster_file, sep='\t', index_col=0 if 'ClusterID' in pd.read_csv(args.cluster_file, sep='\t', nrows=0).columns else None)
    except FileNotFoundError:
        print(f"[ERROR] Cluster file not found: {args.cluster_file}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"[WARNING] Cluster file {args.cluster_file} is empty. Creating empty output.")
        pd.DataFrame().to_csv(args.output_prefix + ".tsv", sep='\t')
        sys.exit(0)

    allele_cols = [c for c in df.columns if c.startswith('Chr') or 'Allele_' in c or c == 'Unplaced_or_Translocated']

    col_suffixes = {}
    for col in allele_cols:
        match = re.search(r'([a-zA-Z]+)$', col)
        if match:
            col_suffixes[col] = match.group(1).upper()

    all_clustered_genes = set()
    for _, row in df.iterrows():
        all_clustered_genes.update(extract_genes_from_row(row, allele_cols))

    verified_rows = []
    rejected_rows = []
    total_genes_checked = 0
    total_genes_mapped_to_og = 0
    total_genes_rescued = 0  # rescued-gene counter

    print(f"[INFO] Starting rigorous orthogonal verification & native rescue...")
    print(f"[INFO] Using {len(allele_cols)} structural columns and a minimum OG support ratio of {args.min_support_ratio}")

    for cluster_id, row in df.iterrows():
        genes = extract_genes_from_row(row, allele_cols)

        if not genes:
            rejected_rows.append(row.to_dict())
            continue

        total_genes_checked += len(genes)

        found_ogs = [og_map.get(g) for g in genes if g in og_map]
        total_genes_mapped_to_og += len(found_ogs)

        if found_ogs:
            main_og = max(set(found_ogs), key=found_ogs.count)
            support_ratio = found_ogs.count(main_og) / len(genes)

            if support_ratio >= args.min_support_ratio:
                row_dict = row.to_dict()
                if df.index.name == 'ClusterID' or cluster_id:
                    row_dict['ClusterID'] = cluster_id
                row_dict['Main_OG'] = main_og
                row_dict['OG_Support_Ratio'] = round(support_ratio, 3)

                missing_genes = og_to_genes[main_og] - all_clustered_genes
                if missing_genes:
                    for mg in missing_genes:
                        target_sub = extract_subgenome_from_id(mg)
                        placed = False

                        if target_sub:
                            for col_name, suffix in col_suffixes.items():
                                if suffix == target_sub:
                                    existing = str(row_dict.get(col_name, 'NA'))
                                    if existing in ('NA', 'nan', ''):
                                        row_dict[col_name] = mg
                                    else:
                                        genes_set = set(existing.split(','))
                                        genes_set.add(mg)
                                        row_dict[col_name] = ",".join(sorted(list(genes_set)))
                                    placed = True
                                    break

                        if not placed:
                            existing_unplaced = str(row_dict.get('Unplaced_or_Translocated', 'NA'))
                            if existing_unplaced in ('NA', 'nan', ''):
                                row_dict['Unplaced_or_Translocated'] = mg
                            else:
                                genes_set = set(existing_unplaced.split(','))
                                genes_set.add(mg)
                                row_dict['Unplaced_or_Translocated'] = ",".join(sorted(list(genes_set)))

                        all_clustered_genes.add(mg)
                        total_genes_rescued += 1

                verified_rows.append(row_dict)
                continue

        row_dict = row.to_dict()
        if df.index.name == 'ClusterID' or cluster_id:
            row_dict['ClusterID'] = cluster_id
        rejected_rows.append(row_dict)

    if total_genes_checked > 0:
        hit_rate = total_genes_mapped_to_og / total_genes_checked
        print(f"[INFO] Global Gene ID Hit Rate against OrthoFinder: {hit_rate:.1%}")
        if hit_rate < 0.1:
            print("[WARNING] Extremely low ID hit rate! The gene IDs in your cluster file do NOT match the IDs in Orthogroups.tsv.")
            print("[WARNING] Please check your FASTA header preprocessing steps.")

    if verified_rows:
        v_df = pd.DataFrame(verified_rows)
        cols = v_df.columns.tolist()
        if 'ClusterID' in cols:
            cols.insert(0, cols.pop(cols.index('ClusterID')))
            v_df = v_df[cols]
        v_df.to_csv(args.output_prefix + ".tsv", sep='\t', index=False)
    else:
        print("[FATAL] 0 clusters passed verification!")
        cols_template = ['ClusterID'] + allele_cols + ['Main_OG', 'OG_Support_Ratio']
        if 'Unplaced_or_Translocated' not in cols_template: cols_template.append('Unplaced_or_Translocated')
        pd.DataFrame(columns=cols_template).to_csv(args.output_prefix + ".tsv", sep='\t', index=False)

    if rejected_rows:
        r_df = pd.DataFrame(rejected_rows)
        cols = r_df.columns.tolist()
        if 'ClusterID' in cols:
            cols.insert(0, cols.pop(cols.index('ClusterID')))
            r_df = r_df[cols]
        r_df.to_csv(args.output_prefix + "_rejected.tsv", sep='\t', index=False)

    print("-" * 60)
    print(f"Verification & Rescue Results:")
    print(f"  Verified Clusters : {len(verified_rows)}")
    print(f"  Rejected Clusters : {len(rejected_rows)}")
    print(f"  Genes Rescued     : {total_genes_rescued} (Smartly placed into subgenomes)")
    print("-" * 60)

if __name__ == '__main__':
    main()
