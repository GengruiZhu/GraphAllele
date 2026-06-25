#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import gzip
import shutil
from collections import defaultdict
from Bio import SeqIO
import subprocess
import re

def safe_extract_id(raw_id):
    """
    Safely strips database prefixes (gene:, transcript:) and
    trailing transcript suffixes (e.g., .1, .2, .t1, .T01),
    preserving the core biological ID and internal dots.
    """
    cleaned = str(raw_id).replace('gene:', '').replace('gene-', '').replace('transcript:', '').strip()
    return re.sub(r'\.[tT]?\d+$', '', cleaned)

def is_target_chrom(chrom_name, target_num):
    """
    Smart Regex Matcher: Dynamically identifies subgenomes without hardcoding A-N.
    Matches Chr1A, Chr01At, Chr1_D1 but safely rejects Chr10A when looking for 1.
    """
    # Regex breakdown:
    # ^Chr0?      -> Starts with 'Chr' or 'Chr0'
    # (\d+)       -> Captures the chromosome integer
    # (.*)$       -> Captures the subgenome suffix
    match = re.match(r'^Chr0?(\d+)(.*)$', chrom_name, re.IGNORECASE)
    if match:
        num = int(match.group(1))
        suffix = match.group(2)
        # Ensure the suffix doesn't start with a number to prevent matching Chr10 to Chr1
        if num == target_num and (not suffix or not suffix[0].isdigit()):
            return True
    return False

def decompress_if_needed(path):
    if path.endswith(".gz"):
        out_path = path[:-3]
        if not os.path.exists(out_path):
            print(f"[INFO] Decompressing {path}...")
            with gzip.open(path, "rt") as fin, open(out_path, "w", encoding='utf-8') as fout:
                shutil.copyfileobj(fin, fout)
        return out_path
    return path

def split_gff_by_chromosome(gff_file, out_dir, target_num=None):
    chrom_lines = defaultdict(list)
    headers = []

    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("##"):
                headers.append(line)
                continue
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) >= 9:
                chrom = parts[0]
                if target_num is None or is_target_chrom(chrom, target_num):
                    chrom_lines[chrom].append(line)

    os.makedirs(out_dir, exist_ok=True)
    for chrom, lines in chrom_lines.items():
        with open(os.path.join(out_dir, f"{chrom}.gff"), "w", encoding='utf-8') as fout:
            fout.writelines(headers)  # Preserve GFF3 headers for strict parsers
            fout.writelines(lines)

    return list(chrom_lines.keys())

def split_fasta_by_chromosome(fasta_file, out_dir, target_num=None):
    os.makedirs(out_dir, exist_ok=True)
    chroms = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom_name = record.id
        if target_num is None or is_target_chrom(chrom_name, target_num):
            chroms.append(chrom_name)
            out_path = os.path.join(out_dir, f"{chrom_name}.fasta")
            with open(out_path, "w") as fout:
                SeqIO.write(record, fout, "fasta")

    return chroms

def clean_fasta_headers(fasta_path):
    """
    Surgical formatting of FASTA headers to guarantee they match the BED file IDs.
    """
    if not os.path.exists(fasta_path): return
    records = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        # Utilize the global safe extractor
        rec.id = safe_extract_id(rec.id)
        rec.description = ""
        records.append(rec)
    SeqIO.write(records, fasta_path, "fasta")

def extract_cds_with_gffread(gff_dir, fasta_dir, cds_dir, log_path, gffread_path="gffread", extract_protein=False):
    os.makedirs(cds_dir, exist_ok=True)
    with open(log_path, "a") as log:
        for fname in os.listdir(gff_dir):
            if not fname.endswith(".gff"): continue
            chrom = fname.replace(".gff", "")
            gff_path = os.path.join(gff_dir, fname)
            fa_path = os.path.join(fasta_dir, f"{chrom}.fasta")

            if not os.path.exists(fa_path): continue

            out_cds = os.path.join(cds_dir, f"{chrom}.cds")
            cmd = [gffread_path, gff_path, "-g", fa_path, "-x", out_cds]

            out_pep = None
            if extract_protein:
                out_pep = os.path.join(cds_dir, f"{chrom}.pep")
                cmd += ["-y", out_pep]

            try:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                # Apply surgical ID trimming immediately after extraction
                clean_fasta_headers(out_cds)
                if out_pep:
                    clean_fasta_headers(out_pep)

                log.write(f"[SUCCESS] {chrom} CDS/PEP extracted and IDs unified.\n")
            except Exception as e:
                log.write(f"[ERROR] {chrom} gffread failed: {str(e)}\n")

def convert_gff_to_bed(gff_dir, bed_dir, log_path, feature_type="gene"):
    """
    Generates the highly sensitive .bed files required by JCVI/MCScanX.
    Applies the same safe_extract_id to ensure a 100% match with the FASTA headers.
    """
    os.makedirs(bed_dir, exist_ok=True)
    with open(log_path, "a") as log:
        for fname in os.listdir(gff_dir):
            if not fname.endswith(".gff"): continue
            chrom = fname.replace(".gff", "")
            gff_path = os.path.join(gff_dir, fname)
            bed_path = os.path.join(bed_dir, f"{chrom}.bed")

            try:
                with open(gff_path, 'r', encoding='utf-8') as f, open(bed_path, "w", encoding='utf-8') as fout:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"): continue

                        parts = line.split("\t")
                        if len(parts) < 9: continue

                        if parts[2] == feature_type:
                            # BED is 0-based index
                            start = str(int(parts[3]) - 1)
                            end = parts[4]

                            gene_id = "Unknown"
                            attributes = parts[8]

                            for field in attributes.split(';'):
                                if field.startswith('ID='):
                                    gene_id = field.split('=', 1)[1]
                                    break
                                elif field.startswith('Name='):
                                    gene_id = field.split('=', 1)[1]

                            # Critical Step: Unify BED ID format with FASTA ID format
                            gene_id = safe_extract_id(gene_id)

                            fout.write(f"{parts[0]}\t{start}\t{end}\t{gene_id}\n")

                log.write(f"[SUCCESS] {chrom} BED generated safely.\n")
            except Exception as e:
                log.write(f"[ERROR] {chrom} BED generation failed: {str(e)}\n")

def main():
    parser = argparse.ArgumentParser(description="Split genome and strictly standardize formats for JCVI synteny pipeline.")
    parser.add_argument('--gff', required=True, help="Input whole genome GFF3")
    parser.add_argument('--fasta', required=True, help="Input whole genome FASTA")
    parser.add_argument('--outdir', default='split_output', help="Working directory")
    parser.add_argument('--gffread', default='gffread', help="Path to gffread tool")
    parser.add_argument('--chr', type=int, help="Target chromosome number to extract")
    parser.add_argument('--protein', '-y', action='store_true', help="Extract protein sequences (.pep)")
    parser.add_argument('--feature', default='gene', help="GFF feature type to build BED coordinates")
    args = parser.parse_args()

    gff_file = decompress_if_needed(args.gff)
    fasta_file = decompress_if_needed(args.fasta)

    gff_dir = os.path.join(args.outdir, "gff")
    fasta_dir = os.path.join(args.outdir, "fasta")
    cds_dir = os.path.join(args.outdir, "cds")
    bed_dir = os.path.join(args.outdir, "bed")
    log_path = os.path.join(args.outdir, "prepare_jcvi.log")

    os.makedirs(args.outdir, exist_ok=True)

    # Wipe old log
    if os.path.exists(log_path):
        os.remove(log_path)

    print(f"[INFO] [Phase 1] Extracting dynamic subgenomes for Chromosome Group {args.chr} from GFF...")
    split_gff_by_chromosome(gff_file, gff_dir, args.chr)

    print(f"[INFO] [Phase 2] Extracting matching FASTA sequences...")
    split_fasta_by_chromosome(fasta_file, fasta_dir, args.chr)

    print("[INFO] [Phase 3] Invoking gffread and structurally aligning IDs for JCVI...")
    extract_cds_with_gffread(gff_dir, fasta_dir, cds_dir, log_path, gffread_path=args.gffread, extract_protein=args.protein)

    print("[INFO] [Phase 4] Building highly-sensitive 0-based BED coordinate files...")
    convert_gff_to_bed(gff_dir, bed_dir, log_path, feature_type=args.feature)

    print(f"[SUCCESS] All inputs cleanly prepared and synchronized in: {args.outdir}")

if __name__ == '__main__':
    main()
