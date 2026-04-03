# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import argparse
import sys

def copy_files(src_dir, dst_dir, extensions):
    """
    Moving JCVI
    """
    os.makedirs(dst_dir, exist_ok=True)
    src_dir = os.path.abspath(src_dir)
    dst_dir = os.path.abspath(dst_dir)
    
    count = 0
    for fname in os.listdir(src_dir):
        if any(fname.endswith(ext) for ext in extensions):
            src = os.path.join(src_dir, fname)
            dst = os.path.join(dst_dir, fname)
            
            if os.path.exists(dst): os.remove(dst)
            os.symlink(src, dst)
            count += 1
    print(f"[INFO] Linked {count} files ({extensions}) to {dst_dir}")

def run_shell_script(script_path, work_dir):
    """
    
    """
    script_abs_path = os.path.abspath(script_path)
    print(f"[INFO] Entering working directory: {work_dir}")
    print(f"[INFO] Preparing to run JCVI script: {script_abs_path}")
    
    
    result = subprocess.run(["bash", script_abs_path], cwd=work_dir)
    if result.returncode != 0:
        raise Exception(f"Shell script execution failed with exit code: {result.returncode}")

def collect_anchors(input_dir, output_dir):
    """
    Harvesting alignment results (.anchors)
    """
    os.makedirs(output_dir, exist_ok=True)
    count = 0
    for root, _, files in os.walk(input_dir):
        for f in files:
            
            if f.endswith(".anchors") and not f.endswith(".lifted.anchors"):
                src = os.path.join(root, f)
                dst = os.path.join(output_dir, f)
                shutil.copy2(src, dst)
                count += 1
    print(f"[INFO] Successfully harvested {count} .anchors files to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="JCVI Pipeline Orchestrator for PolyAlleler")
    parser.add_argument("--cds_dir", required=True, help="Directory containing .cds files")
    parser.add_argument("--bed_dir", required=True, help="Directory containing .bed files")
    parser.add_argument("--jcvi_input", default="jcvi_work", help="Working dir for JCVI process")
    parser.add_argument("--anchors_dir", default="anchors_output", help="Final anchors output dir")
    parser.add_argument("--sh_script", required=True, help="Path to run_ortholog.sh")
    args = parser.parse_args()

    copy_files(args.cds_dir, args.jcvi_input, [".cds"])
    copy_files(args.bed_dir, args.jcvi_input, [".bed"])

    try:
        run_shell_script(args.sh_script, args.jcvi_input)
    except Exception as e:
        print(f"[ERROR] JCVI crashed during execution: {e}")
        sys.exit(1)

    collect_anchors(args.jcvi_input, args.anchors_dir)
    print("[SUCCESS] JCVI pipeline tasks completed successfully.")

if __name__ == "__main__":
    main()
