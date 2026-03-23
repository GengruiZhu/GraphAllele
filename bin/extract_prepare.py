# -*- coding: utf-8 -*-
import os
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description="从全基因组提取特定染色体组的 GFF 和 CDS")
    parser.add_argument("--gff", required=True, help="原始全基因组 GFF3 路径")
    parser.add_argument("--fasta", required=True, help="原始全基因组 FASTA 路径")
    parser.add_argument("--chr", required=True, help="目标染色体编号 (数字, 如 1)")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--gffread", default="gffread", help="gffread 路径")
    
    args = parser.parse_args()

    
    target_chr_prefix = "Chr%s" % args.chr
    output_prefix = os.path.join(args.outdir, "Chr%sA_to_Chr%sN" % (args.chr, args.chr))
    os.makedirs(args.outdir, exist_ok=True)

    
    target_chroms = [f"{target_chr_prefix}{c}" for c in "ABCDEFGHIJKLMN"]

    
    filtered_gff = f"{output_prefix}.gff3"
    print(f"[i] 正在从 {args.gff} 提取同源组 Chr{args.chr}...")
    
    with open(args.gff, "r") as fin, open(filtered_gff, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9 and parts[0] in target_chroms:
                fout.write(line + "\n")

    print(f"GFF 提取完成：{filtered_gff}")

    
    cds_output = f"{output_prefix}.cds.fasta"
    cmd = [
        args.gffread,
        filtered_gff,
        "-g", args.fasta,
        "-x", cds_output
    ]

    print("[i] gffread 提取 CDS...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"CDS 提取成功：{cds_output}")
    else:
        print("[!] gffread 运行失败：")
        print(result.stderr)

if __name__ == "__main__":
    main()
