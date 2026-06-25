#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
import argparse
import sys

def link_or_copy(src, dst):
    """
    Safely create a symlink. If the OS or file system rejects symlinks,
    fallback to a physical copy.
    """
    if os.path.exists(dst):
        os.remove(dst)
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)

def prepare_jcvi_inputs(src_dir, dst_dir, target_type):
    """
    Safely transfer files to the JCVI sandbox.
    Automatically standardizes extensions because JCVI strictly requires
    files to end with exactly '.cds' or '.bed'.
    """
    os.makedirs(dst_dir, exist_ok=True)
    src_dir = os.path.abspath(src_dir)
    dst_dir = os.path.abspath(dst_dir)

    count = 0
    for fname in os.listdir(src_dir):
        # Handle cases where upstream outputs .cds.fasta or .bed
        is_match = False
        if target_type == "cds" and ("cds" in fname.lower() or fname.endswith(".fa") or fname.endswith(".fasta")):
            is_match = True
            # Force the destination name to end strictly with .cds
            base_name = fname.split('.cds')[0].split('.fa')[0]
            dst_name = f"{base_name}.cds"

        elif target_type == "bed" and "bed" in fname.lower():
            is_match = True
            # Force the destination name to end strictly with .bed
            base_name = fname.split('.bed')[0]
            dst_name = f"{base_name}.bed"

        if is_match:
            src = os.path.join(src_dir, fname)
            dst = os.path.join(dst_dir, dst_name)
            link_or_copy(src, dst)
            count += 1

    print(f"[INFO] Synchronized {count} {target_type.upper()} files into JCVI sandbox: {dst_dir}")
    return count

def run_shell_script(script_path, work_dir):
    """
    Executes the bash script containing JCVI pairwise commands.
    Captures output to a log file to keep the terminal clean and aid debugging.
    """
    script_abs_path = os.path.abspath(script_path)
    log_file = os.path.join(work_dir, "jcvi_execution.log")

    print(f"[INFO] Entering working directory: {work_dir}")
    print(f"[INFO] Executing JCVI bash script: {script_abs_path}")
    print(f"[INFO] Logging JCVI standard output to: {log_file}")

    with open(log_file, "w") as f_log:
        result = subprocess.run(["bash", script_abs_path], cwd=work_dir, stdout=f_log, stderr=f_log)

    if result.returncode != 0:
        print(f"[ERROR] JCVI Shell script failed with exit code: {result.returncode}")
        print(f"[ERROR] Please check the log file for details: {log_file}")
        raise Exception("JCVI pipeline execution failed.")

def collect_anchors(input_dir, output_dir):
    """
    Harvest the final .anchors files, ignoring noisy intermediate files.
    """
    os.makedirs(output_dir, exist_ok=True)
    count = 0
    for root, _, files in os.walk(input_dir):
        for f in files:
            if f.endswith(".anchors") and not f.endswith(".lifted.anchors"):
                src = os.path.join(root, f)
                dst = os.path.join(output_dir, f)
                # Only copy if source and destination are different paths
                if os.path.abspath(src) != os.path.abspath(dst):
                    shutil.copy2(src, dst)
                count += 1

    print(f"[INFO] Successfully harvested {count} .anchors files into {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="JCVI Pipeline Orchestrator for GraphAllele")
    parser.add_argument("--cds_dir", required=True, help="Directory containing .cds sequences")
    parser.add_argument("--bed_dir", required=True, help="Directory containing .bed coordinates")
    parser.add_argument("--jcvi_input", required=True, help="Working directory for JCVI process")
    # Anchors dir is now optional. If omitted, it defaults to the jcvi_input directory
    # to perfectly match the expectations of GraphAllele.py Step 5.
    parser.add_argument("--anchors_dir", help="Final anchors output dir (defaults to jcvi_input)")
    parser.add_argument("--sh_script", required=True, help="Path to run_ortholog.sh")
    args = parser.parse_args()

    # Set default output directory to prevent file collision in the root folder
    final_anchors_dir = args.anchors_dir if args.anchors_dir else args.jcvi_input

    cds_count = prepare_jcvi_inputs(args.cds_dir, args.jcvi_input, "cds")
    bed_count = prepare_jcvi_inputs(args.bed_dir, args.jcvi_input, "bed")

    if cds_count == 0 or bed_count == 0:
        print("[ERROR] Input directories are empty or missing correctly formatted files.")
        sys.exit(1)

    try:
        run_shell_script(args.sh_script, args.jcvi_input)
    except Exception as e:
        print(f"[FATAL] JCVI process interrupted: {e}")
        sys.exit(1)

    collect_anchors(args.jcvi_input, final_anchors_dir)
    print("[SUCCESS] JCVI synteny identification completed successfully.")

if __name__ == "__main__":
    main()
