#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import signal
import subprocess
import shutil
from Bio import SeqIO
from Bio.Seq import Seq

def run_intra_group_orthofinder(pep_files_list, outdir, gid, threads):
    of_dir = os.path.join(outdir, "01.5.OrthoFinder_Intra")
    fasta_dir = os.path.join(of_dir, "input_fasta")
    os.makedirs(fasta_dir, exist_ok=True)

    # Check for breakpoint resume
    results_base = os.path.join(fasta_dir, "OrthoFinder")
    if os.path.exists(results_base):
        subdirs = [os.path.join(results_base, d) for d in os.listdir(results_base) if d.startswith("Results_")]
        if subdirs:
            latest = max(subdirs, key=lambda x: os.path.getctime(x))
            og_file_existing = os.path.join(latest, "Orthogroups", "Orthogroups.tsv")
            if os.path.exists(og_file_existing) and os.path.getsize(og_file_existing) > 0:
                print(f"[INFO] OrthoFinder results already exist for {gid}, skipping.")
                return og_file_existing

    print(f"[INFO] [{gid}] Cleaning and Preparing FASTA files...")

    copied_count = 0
    for pep_file in pep_files_list:
        if "ref" not in os.path.basename(pep_file).lower() and pep_file.endswith(".pep"):
            new_name = os.path.basename(pep_file).replace(".pep", ".fa")
            dest_file = os.path.join(fasta_dir, new_name)

            cleaned_records = []
            for record in SeqIO.parse(pep_file, "fasta"):
                seq_str = str(record.seq).upper().replace('*', 'X').replace('.', 'X')
                # Safer method to remove trailing 'X'
                seq_str = seq_str.rstrip('X')

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

    print(f"[INFO] [{gid}] Launching OrthoFinder with TACTICAL SNIPER ({copied_count} files, {safe_threads} threads).")

    cmd = ['orthofinder', '-f', fasta_dir, '-t', str(safe_threads), '-og']

    # Use process group ID (setsid) to ensure all child processes (e.g., diamond, mcl) can be killed
    with open(of_log, "w") as f:
        process = subprocess.Popen(cmd, stdout=f, stderr=f, preexec_fn=os.setsid)

    sniped = False
    target_og_file = None

    while process.poll() is None:
        time.sleep(5)

        # Directly check the output directory instead of parsing log text to improve fault tolerance
        if os.path.exists(results_base):
            subdirs = [os.path.join(results_base, d) for d in os.listdir(results_base) if d.startswith("Results_")]
            if subdirs:
                latest_result_dir = max(subdirs, key=lambda x: os.path.getctime(x))
                possible_tsv = os.path.join(latest_result_dir, "Orthogroups", "Orthogroups.tsv")

                marker_file = os.path.join(latest_result_dir, "Orthogroups", "Orthogroups_UnassignedGenes.tsv")

                if os.path.exists(possible_tsv) and os.path.exists(marker_file):
                    if os.path.getsize(possible_tsv) > 0:
                        print(f"\n[TACTICAL SNIPER] [{gid}] Marker file detected. TSV is safely written to disk!")

                        safe_backup_tsv = os.path.join(of_dir, f"{gid}_Orthogroups_sniped.tsv")
                        shutil.copy2(possible_tsv, safe_backup_tsv)
                        print(f"[TACTICAL SNIPER] [{gid}] Data secured at {safe_backup_tsv}.")

                        print(f"[TACTICAL SNIPER] [{gid}] Initiating graceful shutdown (SIGTERM)...")

                        # Step 1
                        os.killpg(os.getpgid(process.pid), signal.SIGTERM)

                        try:
                            # Allow 10 seconds for OrthoFinder and MCL to exit safely
                            process.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            # Step 2: If the process hangs and refuses to exit, trigger the fallback strict kill mechanism
                            print(f"[TACTICAL SNIPER] [{gid}] Process resistant to SIGTERM. Forcing SIGKILL...")
                            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                            # Ensure the process is completely cleared from the process table
                            process.wait()

                        sniped = True
                        target_og_file = safe_backup_tsv
                        break

    # If ended normally without sniping (meaning -og suddenly worked, or ran too fast)
    if not sniped and process.returncode != 0 and process.returncode != -15 and process.returncode != -9:
        # Note: Filter out return codes -15 (SIGTERM) and -9 (SIGKILL) since we intentionally killed it
        print(f"[ERROR] [{gid}] OrthoFinder crashed unexpectedly. Check log: {of_log}")
        sys.exit(1)

    # If no sniping occurred and it ended normally, we need to locate the file again
    if not sniped and target_og_file is None:
        if os.path.exists(results_base):
            subdirs = [os.path.join(results_base, d) for d in os.listdir(results_base) if d.startswith("Results_")]
            if subdirs:
                latest_result_dir = max(subdirs, key=lambda x: os.path.getctime(x))
                target_og_file = os.path.join(latest_result_dir, "Orthogroups", "Orthogroups.tsv")

    # =========================================================================

    if target_og_file and os.path.exists(target_og_file) and os.path.getsize(target_og_file) > 0:
        print(f"[SUCCESS] [{gid}] Matrix Locked successfully!")
        return target_og_file
    else:
        print(f"[ERROR] [{gid}] Sniper logic failed, file missing or corrupted.")
        sys.exit(1)
