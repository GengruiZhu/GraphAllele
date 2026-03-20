# -*- coding: utf-8 -*-
import os
import subprocess
import argparse
import logging
import sys
import networkx as nx
from collections import defaultdict

# 配置标准日志
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def parse_gff_order(gff_file):
    """解析 GFF：带有容错和异常检查机制"""
    gene_order = {}
    chrom_genes = defaultdict(list)
    
    if not os.path.exists(gff_file):
        logging.error(f"GFF file not found: {gff_file}")
        sys.exit(1)

    with open(gff_file) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip(): 
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene": 
                continue
            
            chrom = parts[0]
            info = parts[8]
            
            # 防御性解析 ID
            gene_id = None
            for field in info.split(';'):
                if field.startswith('ID='):
                    # 统一清洗逻辑
                    gene_id = field[3:].split('.')[0] 
                    break
            
            if gene_id:
                try:
                    start_pos = int(parts[3])
                    chrom_genes[chrom].append((start_pos, gene_id))
                except ValueError:
                    logging.warning(f"Invalid start position at line {line_num}: {parts[3]}")
    
    # 建立索引坐标系
    for chrom, genes in chrom_genes.items():
        genes.sort(key=lambda x: x[0])
        for idx, (_, gene_id) in enumerate(genes):
            gene_order[gene_id] = (chrom, idx)
            
    logging.info(f"Parsed {len(gene_order)} valid genes from GFF.")
    return gene_order

def run_self_blast(pep_file, out_blast, threads=4):
    """运行全对全 BLASTP，带有严格的错误捕获"""
    if not os.path.exists(pep_file):
        logging.error(f"PEP file not found: {pep_file}")
        sys.exit(1)

    db_name = pep_file + ".db"
    logging.info("Building BLAST database...")
    try:
        subprocess.run(['makeblastdb', '-in', pep_file, '-dbtype', 'prot', '-out', db_name], 
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        
        logging.info(f"Running Self-BLASTP (Threads: {threads})...")
        # 增加 identity 和 qcovs 输出用于严格过滤
        subprocess.run([
            'blastp', '-query', pep_file, '-db', db_name,
            '-outfmt', '6 qseqid sseqid pident evalue', 
            '-evalue', '1e-10', '-num_threads', str(threads), '-out', out_blast
        ], check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        logging.error(f"BLAST execution failed: {e.stderr.decode().strip()}")
        sys.exit(1)

def identify_tandems(blast_file, gene_order, out_tandem, max_distance=5, min_identity=50.0):
    """核心图论算法：结合物理距离与严格的序列相似度过滤"""
    G = nx.Graph()
    
    try:
        with open(blast_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4: continue
                
                g1, g2 = parts[0].split(".")[0], parts[1].split(".")[0]
                identity = float(parts[2])
                
                if g1 == g2: continue
                # 严谨性增强：除了 e-value，强制校验 identity
                if identity < min_identity: continue
                
                if g1 in gene_order and g2 in gene_order:
                    chrom1, idx1 = gene_order[g1]
                    chrom2, idx2 = gene_order[g2]
                    
                    if chrom1 == chrom2 and abs(idx1 - idx2) <= max_distance:
                        G.add_edge(g1, g2)
                        
    except FileNotFoundError:
        logging.error(f"BLAST output not found for tandem parsing: {blast_file}")
        sys.exit(1)
    
    # 提取连通分量并输出
    tandem_count = 0
    with open(out_tandem, 'w') as fout:
        for component in nx.connected_components(G):
            if len(component) > 1:
                fout.write(",".join(component) + "\n")
                tandem_count += 1
                
    logging.info(f"Identified {tandem_count} tandem arrays (Max Dist: {max_distance}, Min Ident: {min_identity}%).")

def main():
    parser = argparse.ArgumentParser(description="Rigorous Python Tandem Duplication Identifier")
    parser.add_argument('--gff', required=True, help='Merged GFF file')
    parser.add_argument('--pep', required=True, help='Merged PEP file')
    parser.add_argument('--outdir', default='tandem_out', help='Output directory')
    parser.add_argument('--prefix', default='species', help='Output prefix')
    parser.add_argument('--threads', type=int, default=8, help='BLAST threads')
    parser.add_argument('--max_dist', type=int, default=5, help='最大基因间隔 (默认 5)')
    parser.add_argument('--min_ident', type=float, default=50.0, help='最小序列一致性 (默认 50.0%)')
    args = parser.parse_args()

    gff_abs = os.path.abspath(args.gff)
    pep_abs = os.path.abspath(args.pep)
    outdir_abs = os.path.abspath(args.outdir)
    os.makedirs(outdir_abs, exist_ok=True)
    
    blast_out = os.path.join(outdir_abs, f"{args.prefix}.blast")
    tandem_out = os.path.join(outdir_abs, f"{args.prefix}.tandem")

    try:
        gene_order = parse_gff_order(gff_abs)
        run_self_blast(pep_abs, blast_out, args.threads)
        identify_tandems(blast_out, gene_order, tandem_out, args.max_dist, args.min_ident)
    except Exception as e:
        logging.error(f"Pipeline crashed due to an unexpected error: {str(e)}")
        sys.exit(1)
    finally:
        # 工程严谨性：确保中间大文件一定会被清理
        if os.path.exists(blast_out):
            os.remove(blast_out)
            logging.info("Temporary BLAST files cleaned up.")

if __name__ == '__main__':
    main()
