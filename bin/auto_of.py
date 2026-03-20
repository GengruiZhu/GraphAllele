# -*- coding: utf-8 -*-
import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd, cwd=None):
    """
    General function to execute system commands.
    Catches exceptions and exits to prevent cascading errors downstream.
    """
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Auto OrthoFinder step failed: {' '.join(cmd)}")
        sys.exit(1)

def clean_protein_fasta(raw_pep, clean_pep):
    """
    Protein sequence sanitization engine:
    1. Strip trailing stop codons (*)
    2. Replace internal premature stop codons or invalid placeholders (* or .) with unknown amino acid (X)
    3. Filter out exceedingly short invalid sequences (length < 10 aa)
    """
    print(f"[INFO] Sanitizing protein sequences: {raw_pep}")
    valid_records = []
    
    with open(raw_pep, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Extract sequence and convert to uppercase string
            seq_str = str(record.seq).upper()
            
            # Action 1: Strip trailing asterisk (stop codon)
            if seq_str.endswith('*'):
                seq_str = seq_str[:-1]
            
            # Action 2: Replace internal asterisks or dots with the valid unknown amino acid symbol 'X'
            seq_str = seq_str.replace('*', 'X').replace('.', 'X')
            
            # Action 3: Length fail-safe filtering (discard garbage sequences < 10 aa)
            if len(seq_str) >= 10:
                clean_record = SeqRecord(Seq(seq_str), id=record.id, description="")
                valid_records.append(clean_record)
    
    # Write out the sanitized sequences
    with open(clean_pep, "w") as out_handle:
        SeqIO.write(valid_records, out_handle, "fasta")
    
    print(f"[SUCCESS] Sequence sanitization complete! Retained {len(valid_records)} high-quality protein sequences.")

def run_global_orthofinder(gff, fasta, outdir, threads, start_chr, end_chr):
    """
    Automated OrthoFinder scheduling module.
    Note: Retained as a general bioinformatics module. However, for high-ploidy 
    complex genomes like sugarcane, it is strongly recommended to skip this 
    and run OrthoFinder independently.
    """
    print("\n" + "="*60)
    print("[WARNING] Triggering fully automated OrthoFinder clustering engine!")
    print("[WARNING] Strict notice: For high-ploidy large genomes, this step may take days to weeks.")
    print("[WARNING] Highly recommended to run OrthoFinder independently and mount the results using the -og parameter in the main program!")
    print("="*60 + "\n")
    
    # 1. Create an independent working sandbox directory
    of_dir = os.path.join(outdir, "00.OrthoFinder_Auto")
    os.makedirs(of_dir, exist_ok=True)
    raw_pep = os.path.join(of_dir, "WholeGenome_raw.pep")
    clean_pep = os.path.join(of_dir, "WholeGenome_clean.pep")
    
    # 2. Extract whole-genome protein sequences (using gffread)
    if not os.path.exists(raw_pep):
        print(f"[INFO] Using gffread to translate whole-genome protein sequences from GFF and FASTA...")
        run_cmd(['gffread', gff, '-g', fasta, '-y', raw_pep])
    else:
        print("[INFO] Raw whole-genome protein sequences already exist, skipping extraction.")
        
    # 3. Sequence Sanitization
    if not os.path.exists(clean_pep):
        clean_protein_fasta(raw_pep, clean_pep)
    else:
        print("[INFO] Sanitized high-quality protein sequences already exist, skipping sanitization step.")
        
    # 4. Prepare standard input environment for OrthoFinder
    fasta_dir = os.path.join(of_dir, "input_fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    target_pep = os.path.join(fasta_dir, "Species.fa")
    
    # Create a symlink using absolute path to avoid physical copying and save disk space
    if not os.path.exists(target_pep):
        os.symlink(os.path.abspath(clean_pep), target_pep)
        
    # 5. Launch OrthoFinder core algorithm engine
    print(f"[INFO] Launching OrthoFinder (Allocated physical threads: {threads})...")
    run_cmd(['orthofinder', '-f', fasta_dir, '-t', str(threads)])
    
    # 6. Dynamically locate and return the absolute path of the ultimate orthogroup file
    results_base = os.path.join(fasta_dir, "OrthoFinder")
    if not os.path.exists(results_base):
        print("[ERROR] OrthoFinder failed to generate the result directory! Possible run interruption or Out-Of-Memory (OOM).")
        sys.exit(1)
        
    # Look for the latest output directory starting with Results_
    subdirs = [os.path.join(results_base, d) for d in os.listdir(results_base) if d.startswith("Results_")]
    if not subdirs:
        print("[ERROR] Cannot find the specific Results_xxx directory!")
        sys.exit(1)
        
    latest_result_dir = max(subdirs, key=os.path.getmtime)
    
    # Version-compatible path sniffer (Handles both OrthoFinder v1.x and v2.x structures)
    og_file_v2 = os.path.join(latest_result_dir, "Orthogroups", "Orthogroups.tsv")
    og_file_v1 = os.path.join(latest_result_dir, "Orthogroups.tsv")
    
    if os.path.exists(og_file_v2):
        og_file = og_file_v2
    elif os.path.exists(og_file_v1):
        og_file = og_file_v1
    else:
        print(f"[ERROR] Failed to find Orthogroups.tsv at the expected paths in {latest_result_dir}!")
        sys.exit(1)
        
    print(f"[SUCCESS] Automated OrthoFinder run completed successfully! Locked onto the panoramic orthogroup file: {og_file}")
    return og_file
