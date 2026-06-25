#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import argparse
import shutil

def main():
    parser = argparse.ArgumentParser(description="Extract chromosome-specific GFF and CDS for polyploids")
    parser.add_argument("--gff", required=True, help="Path to whole-genome GFF3")
    parser.add_argument("--fasta", required=True, help="Path to whole-genome FASTA")
    parser.add_argument("--chr", required=True, help="Target chromosome group number (e.g., 1)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--gffread", default="gffread", help="Path to gffread executable")
    # Added the sub_list argument to synchronize with the master script
    parser.add_argument("--sub_list", default="A,B,C,D,E,F,G,H,I,J,K,L,M,N",
                        help="Comma-separated list of subgenome suffixes (e.g., A,B,D or At,Dt)")

    args = parser.parse_args()

    # Step 1: Pre-flight checks for I/O and dependencies
    if not os.path.exists(args.gff):
        print(f"[ERROR] GFF file not found: {args.gff}")
        sys.exit(1)

    if not os.path.exists(args.fasta):
        print(f"[ERROR] FASTA file not found: {args.fasta}")
        sys.exit(1)

    if shutil.which(args.gffread) is None:
        print(f"[ERROR] gffread executable not found in PATH or specified path: {args.gffread}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    output_prefix = os.path.join(args.outdir, f"Group_Chr{args.chr}")

    # Step 2: Dynamic subgenome parsing based on user input
    subgenomes = [sg.strip() for sg in args.sub_list.split(',') if sg.strip()]
    target_chroms = set([f"Chr{args.chr}{sg}" for sg in subgenomes])

    filtered_gff = f"{output_prefix}.gff3"
    print(f"[INFO] Extracting homologous group Chr{args.chr} from {args.gff}...")
    print(f"[INFO] Target chromosomes identified: {', '.join(target_chroms)}")

    # Step 3: Extract GFF with newline fix
    extracted_count = 0
    with open(args.gff, "r", encoding="utf-8") as fin, open(filtered_gff, "w", encoding="utf-8") as fout:
        for line in fin:
            # Preserve header lines
            if line.startswith("#"):
                fout.write(line)
                continue

            parts = line.split("\t")
            # parts[0] is the seqid (chromosome name)
            if len(parts) >= 9 and parts[0] in target_chroms:
                # FIX: 'line' from iterators already contains the trailing newline (\n)
                # Writing line + "\n" would create empty lines that break strict GFF parsers
                fout.write(line)
                extracted_count += 1

    print(f"[INFO] GFF extraction complete: {filtered_gff} ({extracted_count} records extracted)")

    if extracted_count == 0:
        print("[WARNING] No records found for the target chromosomes. Check your --sub_list or FASTA sequence naming format.")

    # Step 4: Extract CDS using gffread
    cds_output = f"{output_prefix}.cds.fasta"
    cmd = [
        args.gffread,
        filtered_gff,
        "-g", args.fasta,
        "-x", cds_output
    ]

    print("[INFO] Running gffread to extract CDS sequences...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Step 5: Validate gffread execution
    if result.returncode == 0:
        print(f"[SUCCESS] CDS extraction successful: {cds_output}")
    else:
        print("[ERROR] gffread execution failed:")
        print(result.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
