# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import os

def load_og(og_file):
    """加载基因到 OG 的映射"""
    og_map = {}
    if not os.path.exists(og_file):
        print(f"[!] Warning: Orthogroups file {og_file} not found.")
        return og_map
    df = pd.read_csv(og_file, sep='\t')
    for _, row in df.iterrows():
        og_id = row['Orthogroup']
        for genes in row[1:].dropna():
            for g in str(genes).split(','):
                # 清洗逻辑：默认去除点号后面的转录本后缀
                clean_id = g.strip().split('.')[0]
                og_map[clean_id] = og_id
    return og_map

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster_file', required=True)
    parser.add_argument('--orthogroups', required=True)
    parser.add_argument('--output_prefix', required=True)
    args = parser.parse_args()

    og_map = load_og(args.orthogroups)
    df = pd.read_csv(args.cluster_file, sep='\t')
    allele_cols = [c for c in df.columns if c.startswith('Chr')]
    
    verified_rows = []
    rejected_rows = []
    
    # --- 智能诊断探针 ---
    sample_cluster_genes = []
    for _, row in df.head(5).iterrows():
        sample_cluster_genes.extend([str(row[c]).split('.')[0] for c in allele_cols if pd.notna(row[c]) and str(row[c]) != 'NA'])
    
    print("\n[DEBUG] --- 基因 ID 匹配智能诊断 ---")
    print(f"  OG 字典里的样本 ID: {list(og_map.keys())[:5]}")
    print(f"  Cluster 里的样本 ID: {sample_cluster_genes[:5]}")
    print("------------------------------------\n")

    for _, row in df.iterrows():
        genes = [str(row[c]).split('.')[0] for c in allele_cols if pd.notna(row[c]) and str(row[c]) != 'NA']
        found_ogs = [og_map.get(g) for g in genes if g in og_map]
        
        if found_ogs:
            main_og = max(set(found_ogs), key=found_ogs.count)
            support_ratio = found_ogs.count(main_og) / len(genes)
            if support_ratio >= 0.5:
                row['Main_OG'] = main_og
                verified_rows.append(row)
                continue
        rejected_rows.append(row)

    # 防崩溃：即使 0 个验证通过，也要输出一个带表头的空表，防止下游断崖式报错
    if verified_rows:
        pd.DataFrame(verified_rows).to_csv(args.output_prefix + ".tsv", sep='\t', index=False)
    else:
        print("[!] 严重警告: 0 个等位基因簇通过验证！已生成空模板矩阵防止报错。请检查上面的 DEBUG 诊断信息！")
        empty_df = pd.DataFrame(columns=df.columns.tolist() + ['Main_OG'])
        empty_df.to_csv(args.output_prefix + ".tsv", sep='\t', index=False)
        
    if rejected_rows:
        pd.DataFrame(rejected_rows).to_csv(args.output_prefix + "_rejected.tsv", sep='\t', index=False)
    
    print(f"[✓] Verify done. Verified: {len(verified_rows)}, Rejected: {len(rejected_rows)}")

if __name__ == '__main__':
    main()
