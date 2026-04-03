# -*- coding: utf-8 -*-
import os
import subprocess
import argparse

def clean_and_rescue_fasta(fasta_file):
    """
    """
    if not os.path.exists(fasta_file):
        return
        
    temp_file = fasta_file + ".tmp"
    with open(fasta_file, 'r') as fin, open(temp_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                parts = line.strip().split()
                if len(parts) > 1:
                    real_id = parts[1]
                else:
                    real_id = line[1:].strip()
                
                safe_id = real_id.replace("SoZg.", "SoZg_")
                fout.write(f">{safe_id}\n")
            else:
                clean_seq = line.strip().replace('.', '').replace('*', '')
                if clean_seq:
                    fout.write(clean_seq + "\n")
                    
    os.replace(temp_file, fasta_file)
    print(f"  -> [+] Data cleaned: {os.path.basename(fasta_file)}")


def main():
    parser = argparse.ArgumentParser(description="Extract GFF,CDS and PEP ")
    parser.add_argument("--gff", required=True, help="GFF3")
    parser.add_argument("--fasta", required=True, help="FASTA")
    parser.add_argument("--chr", required=True, help="(Number,like 1)")
    parser.add_argument("--outdir", required=True, help="output")
    parser.add_argument("--gffread", default="gffread", help="gffread")
    
    args = parser.parse_args()

    target_chr_prefix = "Chr%s" % args.chr
    output_prefix = os.path.join(args.outdir, "Chr%sA_to_Chr%sN" % (args.chr, args.chr))
    os.makedirs(args.outdir, exist_ok=True)

    target_chroms = [f"{target_chr_prefix}{c}" for c in "ABCDEFGHIJKLMN"]

    filtered_gff = f"{output_prefix}.gff3"
    print(f"[i] Already {args.gff} extracting Chr{args.chr}...")
    
    with open(args.gff, "r") as fin, open(filtered_gff, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9 and parts[0] in target_chroms:
                fout.write(line + "\n")

    print(f"GFF Extracting completely:{filtered_gff}")

    cds_output = f"{output_prefix}.cds.fasta"
    pep_output = f"{output_prefix}.pep"
    
    cmd = [
        args.gffread,
        filtered_gff,
        "-g", args.fasta,
        "-x", cds_output,
        "-y", pep_output
    ]

    print("[i] Beginning gffread extracting CDS and PEP ...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print(f"[√] Successuly")
        clean_and_rescue_fasta(cds_output)
        clean_and_rescue_fasta(pep_output)
        print("[√] All done!")
    else:
        print("[!] gffread 运Failed")
        print(result.stderr)

if __name__ == "__main__":
    main()
