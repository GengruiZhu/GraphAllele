# -*- coding: utf-8 -*-
import os
import sys
import time
import signal
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

def run_cmd(cmd, cwd=None, log_file=None):
    """Log"""
    try:
        if log_file:
            with open(log_file, "w") as f:
                subprocess.run(cmd, check=True, cwd=cwd, stdout=f, stderr=f)
        else:
            subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"\n[FATAL ERROR] External software crashed!")
        print(f"[COMMAND] {' '.join(cmd)}")
        print(f"[CHECK LOG] {log_file if log_file else 'N/A'}")
        sys.exit(1)

def run_intra_group_orthofinder(pep_files_list, outdir, gid, threads):
    of_dir = os.path.join(outdir, "01.5.OrthoFinder_Intra")
    fasta_dir = os.path.join(of_dir, "input_fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    
    
    results_base = os.path.join(fasta_dir, "OrthoFinder")
    if os.path.exists(results_base):
        subdirs = [os.path.join(results_base, d) for d in os.listdir(results_base) if d.startswith("Results_")]
        if subdirs:
            latest = max(subdirs, key=os.path.getmtime)
            og_file_v2 = os.path.join(latest, "Orthogroups", "Orthogroups.tsv")
            if os.path.exists(og_file_v2):
                print(f"[INFO] OrthoFinder results already exist for {gid}, skipping.")
                return og_file_v2

    print(f"[INFO] [{gid}] Cleaning and Preparing FASTA files...")
    
    copied_count = 0
    for pep_file in pep_files_list:
        if "ref" not in os.path.basename(pep_file).lower() and pep_file.endswith(".pep"):
            new_name = os.path.basename(pep_file).replace(".pep", ".fa")
            dest_file = os.path.join(fasta_dir, new_name)
            
            cleaned_records = []
            for record in SeqIO.parse(pep_file, "fasta"):
                seq_str = str(record.seq).upper().replace('*', 'X').replace('.', 'X')
                if seq_str.endswith('X'): seq_str = seq_str[:-1]
                if len(seq_str) >= 10:
                    record.seq = Seq(seq_str)
                    cleaned_records.append(record)
            
            if cleaned_records:
                SeqIO.write(cleaned_records, dest_file, "fasta")
                copied_count += 1
            
    if copied_count < 2:
        print(f"[ERROR] [{gid}] Only {copied_count} species found. Need >= 2.")
        sys.exit(1)

    safe_threads = threads
    of_log = os.path.join(of_dir, f"{gid}_orthofinder_run.log")

    print(f"[INFO] [{gid}] Launching OrthoFinder with VISUAL SNIPER ({copied_count} files, {safe_threads} threads).")
    
    cmd = ['orthofinder', '-f', fasta_dir, '-t', str(safe_threads), '-og']
    
    
    
    
    with open(of_log, "w") as f:
        process = subprocess.Popen(cmd, stdout=f, stderr=f, preexec_fn=os.setsid)
        
    sniped = False
    target_og_file = None

    while process.poll() is None:
        time.sleep(5)  
        
        
        log_ready = False
        if os.path.exists(of_log):
            with open(of_log, "r") as log_f:
                log_content = log_f.read()
                if "Done orthogroups" in log_content or "Starting MSA/Trees" in log_content:
                    log_ready = True

        
        if log_ready:
            if os.path.exists(results_base):
                subdirs = [os.path.join(results_base, d) for d in os.listdir(results_base) if d.startswith("Results_")]
                if subdirs:
                    latest_result_dir = max(subdirs, key=os.path.getmtime)
                    possible_tsv = os.path.join(latest_result_dir, "Orthogroups", "Orthogroups.tsv")
                    
                    
                    if os.path.exists(possible_tsv) and os.path.getsize(possible_tsv) > 0:
                        print(f"\n[SNIPER] [{gid}] Log confirmed AND Target TSV visually confirmed on disk!")
                        print(f"[SNIPER] [{gid}] Waiting 10 seconds for OS to safely close the file handle...")
                        time.sleep(10)
                        
                        print(f"[SNIPER] [{gid}] Assassinating OrthoFinder to prevent combinatorial explosion...")
                        os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                        sniped = True
                        target_og_file = possible_tsv
                        break
                    
    
    if not sniped and process.returncode != 0:
        print(f"[ERROR] [{gid}] OrthoFinder crashed before TSV generation. Check log: {of_log}")
        sys.exit(1)
        
    # =========================================================================
    
    if target_og_file and os.path.exists(target_og_file):
        print(f"[SUCCESS] [{gid}] Matrix Locked successfully by Visual Sniper!")
        return target_og_file
    else:
        print(f"[ERROR] [{gid}] Sniper logic failed, file not found after kill.")
        sys.exit(1)
