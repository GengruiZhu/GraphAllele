# -*- coding: utf-8 -*-
import os
import argparse
import gzip
import shutil
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq  # Added import: for reconstructing cleaned sequence objects
import subprocess

def decompress_if_needed(path):
    if path.endswith(".gz"):
        out_path = path[:-3]
        if not os.path.exists(out_path):
            print(f"[i] Decompressing {path}...")
            with gzip.open(path, "rt") as fin, open(out_path, "w") as fout:
                shutil.copyfileobj(fin, fout)
        return out_path
    return path

def split_gff_by_chromosome(gff_file, out_dir, target_chroms=None):
    chrom_lines = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): 
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9:
                chrom = parts[0]
                if not target_chroms or chrom in target_chroms:
                    chrom_lines[chrom].append(line)
    
    os.makedirs(out_dir, exist_ok=True)
    for chrom, lines in chrom_lines.items():
        with open(os.path.join(out_dir, "%s.gff" % chrom), "w") as fout:
            fout.writelines(lines)
    return list(chrom_lines.items())

def split_fasta_by_chromosome(fasta_file, out_dir, target_chroms=None):
    os.makedirs(out_dir, exist_ok=True)
    chroms = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom_name = record.id
        if not target_chroms or chrom_name in target_chroms:
            chroms.append(chrom_name)
            out_path = os.path.join(out_dir, "%s.fasta" % chrom_name)
            with open(out_path, "w") as fout:
                SeqIO.write(record, fout, "fasta")
    return chroms

def clean_fasta_headers(fasta_path):
    """
    [Ultimate Armor Cleaning Module]
    Fixes ID truncation and cleans illegal characters (*, ., -, #)
    """
    if not os.path.exists(fasta_path): return
    records = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        # --- 1. Rescue ID ---
        # Disguise dangerous 'SoZg.' as 'SoZg_' to prevent accidental truncation
        safe_id = rec.id.replace("SoZg.", "SoZg_")
        # Safely trim real transcript suffixes like .t1
        rec.id = safe_id.split('.')[0]
        rec.description = ""
        
        # --- 2. Clean Sequence ---
        # Forcefully eradicate all abnormal characters that could crash BLAST
        clean_seq_str = str(rec.seq).replace('*', '').replace('.', '').replace('-', '').replace('#', '')
        
        # Extreme case defense: if sequence becomes completely empty after cleaning, discard it
        if len(clean_seq_str) == 0:
            continue
            
        rec.seq = Seq(clean_seq_str)
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
                subprocess.run(cmd, check=True)
                # Execute "plastic surgery" immediately after extraction
                clean_fasta_headers(out_cds)
                if out_pep:
                    clean_fasta_headers(out_pep)
                    
                log.write(f"[✓] {chrom} CDS/PEP extracted and IDs aligned (Cleaned)\n")
            except Exception as e:
                log.write(f"[✗] {chrom} gffread failed: {str(e)}\n")

def convert_gff_to_bed(gff_dir, bed_dir, log_path, feature_type="gene"):
    os.makedirs(bed_dir, exist_ok=True)
    with open(log_path, "a") as log:
        for fname in os.listdir(gff_dir):
            if not fname.endswith(".gff"): continue
            chrom = fname.replace(".gff", "")
            gff_path = os.path.join(gff_dir, fname)
            bed_path = os.path.join(bed_dir, f"{chrom}.bed")
            
            try:
                with open(gff_path) as f, open(bed_path, "w") as fout:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"): continue
                        
                        parts = line.split("\t")
                        if len(parts) < 9: continue
                        
                        if parts[2] == feature_type:
                            start = str(int(parts[3]) - 1)
                            end = parts[4]
                            gene_id = "Unknown"
                            attributes = parts[8]
                            if "ID=" in attributes:
                                raw_id = attributes.split("ID=")[1].split(";")[0]
                            elif "Name=" in attributes:
                                raw_id = attributes.split("Name=")[1].split(";")[0]
                            else:
                                raw_id = gene_id
                                
                            safe_id = raw_id.replace("SoZg.", "SoZg_")
                            final_id = safe_id.split('.')[0]
                                
                            fout.write(f"{parts[0]}\t{start}\t{end}\t{final_id}\n")
                                
                log.write(f"[✓] {chrom} BED written\n")
            except Exception as e:
                log.write(f"[✗] {chrom} BED failed: {str(e)}\n")

def get_chr_list(chr_num):
    letters = "ABCDEFGHIJKLMN"
    list_v1 = ["Chr%d%s" % (chr_num, c) for c in letters]
    list_v2 = ["Chr%02d%s" % (chr_num, c) for c in letters]
    return list(set(list_v1 + list_v2))

def main():
    parser = argparse.ArgumentParser(description="Split GFF/FASTA and extract CDS/PEP/BED per chromosome")
    parser.add_argument('--gff', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--outdir', default='split_output')
    parser.add_argument('--gffread', default='gffread')
    parser.add_argument('--chr', type=int)
    parser.add_argument('--protein', '-y', action='store_true')
    parser.add_argument('--feature', default='gene')
    args = parser.parse_args()

    gff_file = decompress_if_needed(args.gff)
    fasta_file = decompress_if_needed(args.fasta)

    target_chroms = get_chr_list(args.chr) if args.chr else None

    gff_dir = os.path.join(args.outdir, "gff")
    fasta_dir = os.path.join(args.outdir, "fasta")
    cds_dir = os.path.join(args.outdir, "cds")
    bed_dir = os.path.join(args.outdir, "bed")
    log_path = os.path.join(args.outdir, "log.txt")
    os.makedirs(args.outdir, exist_ok=True)

    print("[1] Splitting GFF...")
    split_gff_by_chromosome(gff_file, gff_dir, target_chroms)

    print("[2] Splitting FASTA...")
    split_fasta_by_chromosome(fasta_file, fasta_dir, target_chroms)

    print("[3] Extracting CDS and PEP & Aligning IDs (with Deep Cleaning)...")
    extract_cds_with_gffread(gff_dir, fasta_dir, cds_dir, log_path, gffread_path=args.gffread, extract_protein=args.protein)

    print("[4] Converting to BED...")
    convert_gff_to_bed(gff_dir, bed_dir, log_path, feature_type=args.feature)

    print("[✓] All done. Outputs in: %s" % args.outdir)

if __name__ == '__main__':
    main()
