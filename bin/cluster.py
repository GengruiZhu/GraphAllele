# -*- coding: utf-8 -*-
import os, glob, argparse
import pandas as pd
import networkx as nx
from collections import defaultdict

def parse_gff(gff_file):
    gene_pos = {}
    chrom_order = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene": continue
            
            gene_id = cols[8].split("ID=")[1].split(";")[0].split(".")[0]
            chrom_order[cols[0]].append(gene_id)
            
            gene_pos[gene_id] = (cols[0], len(chrom_order[cols[0]]) - 1)
    return gene_pos

def load_tandem(tandem_file):
    if not tandem_file or not os.path.exists(tandem_file): return set()
    blacklist = set()
    with open(tandem_file) as f:
        for line in f:
            if line.startswith("#"): continue
            blacklist.update([g.strip().split(".")[0] for g in line.split(",")])
    return blacklist

def run_clustering(gff_file, jcvi_dir, tandem_file, output_file, min_chroms, max_dist, sub_list, chr_num):
    gene_pos = parse_gff(gff_file)
    blacklist = load_tandem(tandem_file)
    G = nx.Graph()
    
    
    for fn in glob.glob(os.path.join(jcvi_dir, "*.anchors")):
        with open(fn) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split()
                g1, g2 = [x.split(".")[0] for x in parts[:2]]
                if g1 in blacklist or g2 in blacklist: continue 
                if g1 in gene_pos and g2 in gene_pos:
                    c1, i1 = gene_pos[g1][0], gene_pos[g1][1]
                    c2, i2 = gene_pos[g2][0], gene_pos[g2][1]
                    
                    if abs(i1 - i2) <= max_dist or c1 != c2: G.add_edge(g1, g2)

    
    subs = sub_list.split(',')
    standard_cols = [f"Chr{int(chr_num):02d}{s}" for s in subs]
    
    
    results = []
    for component in nx.connected_components(G):
        
        if len(set(gene_pos[g][0] for g in component)) >= min_chroms:
            row_map = {col: "NA" for col in standard_cols}
            
            for gene in component:
                
                chrom_name = gene_pos[gene][0]
                
                
                for s, col_name in zip(subs, standard_cols):
                    
                    if chrom_name.endswith(s):
                        row_map[col_name] = gene
                        break
                        
            results.append(row_map)

    
    df = pd.DataFrame(results, columns=standard_cols)
    df.to_csv(output_file, sep='\t', index_label="ClusterID")
    print(f"[✓] Standardized cluster done: {len(results)} groups.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', required=True)
    parser.add_argument('--jcvi_dir', required=True)
    parser.add_argument('--tandem', help='Tandem blacklist')
    parser.add_argument('--output', required=True)
    parser.add_argument('--min_chroms', type=int, default=3)
    parser.add_argument('--max_gene_distance', type=int, default=30)
    parser.add_argument('--sub_list', required=True)
    parser.add_argument('--chr_num', required=True)
    args = parser.parse_args()
    run_clustering(args.gff, args.jcvi_dir, args.tandem, args.output, args.min_chroms, args.max_gene_distance, args.sub_list, args.chr_num)
