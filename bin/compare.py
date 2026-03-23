# -*- coding: utf-8 -*-
import pandas as pd
import argparse, subprocess, os
from Bio import SeqIO

def parse_gff_locus(gff_file):
    """"""
    locus_map = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.split("\t")
            if parts[2] == "gene":
                
                try:
                    gene_id = parts[8].split("ID=")[1].split(";")[0]
                    locus = f"{parts[0]}:{parts[3]}-{parts[4]}({parts[6]})"
                    locus_map[gene_id.split('.')[0]] = locus
                except IndexError:
                    continue
    return locus_map

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--allele_file', required=True)
    parser.add_argument('--ref_gff', required=True)
    parser.add_argument('--ref_cds', required=True)
    parser.add_argument('--hap_cds', required=True) 
    parser.add_argument('--output', required=True)
    parser.add_argument('--sub_list', required=True)
    parser.add_argument('--chr_num', required=True)
    parser.add_argument('--identity', type=float, default=80.0) 
    parser.add_argument('--threads', type=int, default=10)
    args = parser.parse_args()

    
    subs = args.sub_list.split(',')
    standard_cols = [f"Chr{int(args.chr_num):02d}{s}" for s in subs]
    
    
    locus_map = parse_gff_locus(args.ref_gff)

    
    db_name = args.ref_cds + ".db"
    
    if not os.path.exists(db_name + ".nin"):
        subprocess.run(['makeblastdb', '-in', args.ref_cds, '-dbtype', 'nucl', '-out', db_name], check=True, stdout=subprocess.DEVNULL)
    
    blast_tmp = args.output + ".blast.tmp"
    
    
    print("[INFO] running tblastn ...")
    subprocess.run(['tblastn', '-query', args.hap_cds, '-db', db_name, '-outfmt', '6 qseqid sseqid pident', 
                    '-num_threads', str(args.threads), '-out', blast_tmp], check=True)

    
    best_hits = {}
    
    if os.path.exists(blast_tmp) and os.path.getsize(blast_tmp) > 0:
        b_df = pd.read_csv(blast_tmp, sep='\t', header=None)
        b_df = b_df[b_df[2] >= args.identity]
        for _, r in b_df.iterrows():
            q_id = str(r[0]).split('.')[0]
            s_id = str(r[1]).split('.')[0]
            if q_id not in best_hits: best_hits[q_id] = s_id
    else:
        print("[!] TBLASTN NALL")

    
    df = pd.read_csv(args.allele_file, sep='\t')
    final_rows = []
    
    for _, row in df.iterrows():
        
        ref_id = "NA"
        ref_locus = "NA"
        
        
        cluster_id = row['ClusterID'] if 'ClusterID' in row else f"Cluster_{len(final_rows)}"
        
        for col in standard_cols:
            if col in row:
                gene = str(row[col]).split('.')[0]
                if gene in best_hits:
                    ref_id = best_hits[gene]
                    ref_locus = locus_map.get(ref_id, "NA")
                    break
        
        
        new_row = [cluster_id, ref_id, ref_locus] + [row.get(col, "NA") for col in standard_cols]
        final_rows.append(new_row)

    
    out_cols = ['ClusterID', 'Ref_Gene', 'Ref_Locus'] + standard_cols
    final_df = pd.DataFrame(final_rows, columns=out_cols)
    final_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"[✓] Final Standardized Table generated: {args.output}")
    
    
    if os.path.exists(blast_tmp): os.remove(blast_tmp)

if __name__ == '__main__':
    main()
