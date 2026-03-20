# -*- coding: utf-8 -*-
import pandas as pd
import argparse, subprocess, os
from Bio import SeqIO

def parse_gff_locus(gff_file):
    """提取参考基因组的坐标信息"""
    locus_map = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.split("\t")
            if parts[2] == "gene":
                # 兼容不同的 GFF 格式
                try:
                    gene_id = parts[8].split("ID=")[1].split(";")[0]
                    locus = f"{parts[0]}:{parts[3]}-{parts[4]}({parts[6]})"
                    locus_map[gene_id.split('.')[0]] = locus
                except IndexError:
                    continue
    return locus_map

def main():
    parser = argparse.ArgumentParser(description="终审校准：对比外部参考基因组")
    parser.add_argument('--allele_file', required=True)
    parser.add_argument('--ref_gff', required=True)
    parser.add_argument('--ref_cds', required=True)
    parser.add_argument('--hap_cds', required=True) # 这里接收到的是 m_pep
    parser.add_argument('--output', required=True)
    parser.add_argument('--sub_list', required=True)
    parser.add_argument('--chr_num', required=True)
    parser.add_argument('--identity', type=float, default=80.0) # 考虑到蛋白比对翻译核酸，阈值稍微放宽
    parser.add_argument('--threads', type=int, default=10)
    args = parser.parse_args()

    # 1. 准备标准列名模板
    subs = args.sub_list.split(',')
    standard_cols = [f"Chr{int(args.chr_num):02d}{s}" for s in subs]
    
    # 2. 提取参考基因组坐标
    locus_map = parse_gff_locus(args.ref_gff)

    # 3. 运行 TBLASTN (Query: 蛋白单倍型, Subject: 外部参考核酸)
    db_name = args.ref_cds + ".db"
    # 如果数据库不存在才建库，节省时间
    if not os.path.exists(db_name + ".nin"):
        subprocess.run(['makeblastdb', '-in', args.ref_cds, '-dbtype', 'nucl', '-out', db_name], check=True, stdout=subprocess.DEVNULL)
    
    blast_tmp = args.output + ".blast.tmp"
    
    # 【核心修复 1】：将 blastn 改为 tblastn
    print("[INFO] 正在运行 tblastn 对齐参考基因组...")
    subprocess.run(['tblastn', '-query', args.hap_cds, '-db', db_name, '-outfmt', '6 qseqid sseqid pident', 
                    '-num_threads', str(args.threads), '-out', blast_tmp], check=True)

    # 4. 读取比对结果：建立 我们的基因 -> 参考基因 的映射
    best_hits = {}
    # 【核心修复 2】：不仅判断文件存在，还要判断文件不为空 (size > 0)
    if os.path.exists(blast_tmp) and os.path.getsize(blast_tmp) > 0:
        b_df = pd.read_csv(blast_tmp, sep='\t', header=None)
        b_df = b_df[b_df[2] >= args.identity]
        for _, r in b_df.iterrows():
            q_id = str(r[0]).split('.')[0]
            s_id = str(r[1]).split('.')[0]
            if q_id not in best_hits: best_hits[q_id] = s_id
    else:
        print("[!] 警告：TBLASTN 结果为空，无法锚定参考基因组！")

    # 5. 读取等位基因表并补全参考信息
    df = pd.read_csv(args.allele_file, sep='\t')
    final_rows = []
    
    for _, row in df.iterrows():
        # 寻找该簇对应的参考基因
        ref_id = "NA"
        ref_locus = "NA"
        
        # 确保 ClusterID 列存在
        cluster_id = row['ClusterID'] if 'ClusterID' in row else f"Cluster_{len(final_rows)}"
        
        for col in standard_cols:
            if col in row:
                gene = str(row[col]).split('.')[0]
                if gene in best_hits:
                    ref_id = best_hits[gene]
                    ref_locus = locus_map.get(ref_id, "NA")
                    break
        
        # 组装最终行
        new_row = [cluster_id, ref_id, ref_locus] + [row.get(col, "NA") for col in standard_cols]
        final_rows.append(new_row)

    # 6. 输出标准矩阵
    out_cols = ['ClusterID', 'Ref_Gene', 'Ref_Locus'] + standard_cols
    final_df = pd.DataFrame(final_rows, columns=out_cols)
    final_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"[✓] Final Standardized Table generated: {args.output}")
    
    # 清理临时文件
    if os.path.exists(blast_tmp): os.remove(blast_tmp)

if __name__ == '__main__':
    main()
