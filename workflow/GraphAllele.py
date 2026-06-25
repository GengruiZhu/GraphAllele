#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GraphAllele - Standardized Allele Matrix Pipeline (Graph-Constrained)

Function: Supports preset subgenome structures, breakpoint resume, dynamic graph
          constraints, and rigorous ID synchronization for downstream workflows.
Notes:    Integrated Smart Subgenome Backfilling and intra-group Auto-OrthoFinder
          with early termination (Tactical Sniper).
"""

import os
import sys
import subprocess
import argparse
import shutil
import datetime
import glob
import re

WORKFLOW_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(WORKFLOW_DIR)
BIN_DIR = os.path.join(PROJECT_ROOT, "bin")

def get_now():
    """Return the current timestamp in a standardized format."""
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def run_cmd(cmd, cwd=None):
    """Execute a shell command strictly with subprocess and error handling."""
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[{get_now()}] [ERROR] Command failed: {' '.join(map(str, cmd))}")
        return False
    except (OSError, TypeError) as e:
        print(f"[{get_now()}] [ERROR] Cannot execute command: {' '.join(map(str, cmd))}")
        print(f"[{get_now()}] [ERROR] {e}")
        return False

def require_file(path, label):
    """Fail early with a clear message if a required input/output file is absent."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        details = [
            f"{label} missing or empty: {path}",
            f"cwd: {os.getcwd()}",
            f"absolute candidate: {os.path.abspath(path) if path else 'NA'}",
        ]
        if path:
            details.append(f"lexists: {os.path.lexists(path)}")
            details.append(f"exists: {os.path.exists(path)}")
            if os.path.islink(path):
                details.append(f"symlink target: {os.readlink(path)}")
                details.append(f"realpath: {os.path.realpath(path)}")
            if os.path.exists(path) and os.path.isfile(path):
                details.append(f"size: {os.path.getsize(path)}")
        raise FileNotFoundError("; ".join(details))

def step_done(marker, required_paths=None, allow_empty_paths=None):
    """Treat a step as complete only when its marker and critical outputs both exist."""
    if not os.path.exists(marker):
        return False
    required_paths = required_paths or []
    allow_empty_paths = set(os.path.abspath(p) for p in (allow_empty_paths or []))
    missing = []
    for p in required_paths:
        if not os.path.exists(p):
            missing.append(p)
        elif os.path.isdir(p) and not os.listdir(p):
            missing.append(p)
        elif os.path.isfile(p) and os.path.getsize(p) == 0 and os.path.abspath(p) not in allow_empty_paths:
            missing.append(p)
    if missing:
        print(f"[{get_now()}] [WARNING] Found {marker}, but required outputs are missing. Rerunning step.")
        for p in missing:
            print(f"[{get_now()}] [WARNING] Missing output: {p}")
        try:
            os.remove(marker)
        except OSError:
            pass
        return False
    return True

def step_done_with_glob(marker, required_paths=None, required_patterns=None):
    """Extended completion check for steps whose outputs are best expressed as globs."""
    if not step_done(marker, required_paths or []):
        return False
    for pattern in required_patterns or []:
        hits = [p for p in glob.glob(pattern) if os.path.isfile(p) and os.path.getsize(p) > 0]
        if not hits:
            print(f"[{get_now()}] [WARNING] Found {marker}, but no output matched: {pattern}. Rerunning step.")
            try:
                os.remove(marker)
            except OSError:
                pass
            return False
    return True

def generate_clean_clusters_auto(work_dir, sub_list_str, out_dense="PolyAlleler_Global_Matrix_Cleaned.tsv", out_kmer="my_clusters.tsv"):
    """
    Core ID Synchronization Module (Smart Backfilling Version):
    1. Pristine matrix mapping valid alleles to their subgenome columns.
    2. Smart Regex Engine: Auto-assigns salvaged genes back to subgenome columns if their ID contains subgenome signatures.
    3. Catch-all 'Salvaged_Genes' column for truly unplaceable contigs.
    """
    print("\n" + "="*60)
    print(f"[{get_now()}] [INFO] Triggering matrix cleaning and Smart Subgenome Backfilling...")

    search_pattern = os.path.join(work_dir, "Group_Chr*", "07.FINAL_ALLELE.tsv")
    file_list = glob.glob(search_pattern)

    if not file_list:
        print(f"[{get_now()}] [WARNING] No 07.FINAL_ALLELE.tsv files found in {work_dir}. Skipping.")
        return

    path_dense = os.path.join(work_dir, out_dense)
    path_kmer = os.path.join(work_dir, out_kmer)

    subgenomes = [sg.strip() for sg in sub_list_str.split(',') if sg.strip()]
    global_idx = 0
    total_salvaged = 0
    auto_backfilled = 0

    null_vals = ('-', '.', 'na', 'nan', '')

    with open(path_dense, 'w', encoding='utf-8') as fdense, open(path_kmer, 'w', encoding='utf-8') as fkmer:

        # Update Header: Add the 'Salvaged_Genes' catch-all column at the very end
        header_cols = ['ClusterID', 'Ref_Gene', 'Ref_Locus'] + [f'Allele_{sg}' for sg in subgenomes] + ['Salvaged_Genes']
        fdense.write('\t'.join(header_cols) + '\n')

        for file_path in sorted(file_list):
            with open(file_path, 'r', encoding='utf-8') as fin:
                lines = fin.readlines()
                if not lines: continue

                raw_header = lines[0].strip('\n').split('\t')
                sg_map = {}
                unplaced_idx = -1

                for i in range(3, len(raw_header)):
                    col_name = raw_header[i].strip()
                    if col_name == "Unplaced_or_Translocated":
                        unplaced_idx = i
                        continue
                    for sg in subgenomes:
                        if col_name.endswith(sg):
                            sg_map[i] = sg
                            break

                for line in lines[1:]:
                    line = line.strip('\n')
                    if not line or line.startswith('#'): continue
                    parts = line.split('\t')
                    if len(parts) < 3: continue

                    ref_g = parts[1].strip()
                    ref_l = parts[2].strip()

                    allele_dict = {sg: 'NA' for sg in subgenomes}
                    flat_alleles = []
                    leftover_salvaged = []

                    # Phase A: Process standard synteny columns
                    for i in range(3, len(parts)):
                        if i == unplaced_idx: continue

                        val = parts[i].strip()
                        if str(val).lower() not in null_vals:
                            genes_in_cell = [v.strip() for v in val.split(',') if v.strip() and str(v.strip()).lower() not in null_vals]
                            if not genes_in_cell: continue

                            if i in sg_map:
                                sg = sg_map[i]
                                allele_dict[sg] = ",".join(genes_in_cell) if allele_dict[sg] == 'NA' else allele_dict[sg] + "," + ",".join(genes_in_cell)

                            flat_alleles.extend(genes_in_cell)

                    # Phase B: Process Salvaged column with SMART BACKFILLING
                    if unplaced_idx != -1 and len(parts) > unplaced_idx:
                        val = parts[unplaced_idx].strip()
                        if str(val).lower() not in null_vals:
                            salvaged_genes = [v.strip() for v in val.split(',') if v.strip() and str(v.strip()).lower() not in null_vals]

                            for gene in salvaged_genes:
                                total_salvaged += 1
                                assigned_sg = None

                                # Regex logic: Look for 'Chr' + any digits + 'Subgenome' NOT followed by another letter.
                                # Matches e.g. <prefix>_Chr1F0001 (F), <prefix>_Chr01E_0001 (E) safely.
                                for sg in subgenomes:
                                    pattern = r'Chr\d+' + re.escape(sg) + r'(?![a-zA-Z])'
                                    if re.search(pattern, gene, re.IGNORECASE):
                                        assigned_sg = sg
                                        break

                                if assigned_sg:
                                    # Backfill into the correct subgenome column
                                    allele_dict[assigned_sg] = gene if allele_dict[assigned_sg] == 'NA' else allele_dict[assigned_sg] + f",{gene}"
                                    auto_backfilled += 1
                                else:
                                    # Truly unplaceable (e.g., Scaffold_123_gene1)
                                    leftover_salvaged.append(gene)

                                flat_alleles.append(gene)

                    # Phase C: Write outputs
                    if flat_alleles:
                        global_id = f"Global_Cluster_{global_idx:06d}"

                        # Serialize leftover salvaged genes
                        salvaged_str = ",".join(leftover_salvaged) if leftover_salvaged else "NA"

                        # Write to dense matrix
                        row = [global_id, ref_g, ref_l] + [allele_dict[sg] for sg in subgenomes] + [salvaged_str]
                        fdense.write('\t'.join(row) + '\n')

                        # Write to flat list (clean, no nan)
                        fkmer.write(f"{global_id}\t{','.join(flat_alleles)}\n")
                        global_idx += 1

    print(f"[{get_now()}] [SUCCESS] Matrix generated! {global_idx} Global Clusters securely formatted.")
    print(f"[{get_now()}] [STAT] Total Salvaged Genes Processed: {total_salvaged}")
    print(f"[{get_now()}] [STAT] Successfully Backfilled to Subgenomes: {auto_backfilled}")
    print(f"[{get_now()}] [STAT] Leftover unplaceable contigs: {total_salvaged - auto_backfilled}")
    print(f"[{get_now()}] [SUCCESS] Pristine Matrix: {path_dense}")
    print(f"[{get_now()}] [SUCCESS] Complete Target List: {path_kmer}")
    print("="*60 + "\n")

def run_group_pipeline(args_tuple):
    chr_num, g_args = args_tuple
    gid = "Group_Chr%02d" % chr_num
    g_dir = os.path.abspath(os.path.join(g_args.outdir, gid))
    os.makedirs(g_dir, exist_ok=True)

    f_abs = os.path.abspath(g_args.fasta)
    g_abs = os.path.abspath(g_args.gff)
    ref_f_abs = os.path.abspath(g_args.ref_cds)
    ref_g_abs = os.path.abspath(g_args.ref_gff)
    og_abs = os.path.abspath(g_args.orthogroups) if g_args.orthogroups else None

    try:
        # Step 1: Prepare Data
        s1_out = os.path.join(g_dir, "01.prepare")
        s1_done = os.path.join(g_dir, ".step1_done")
        s1_required = [os.path.join(s1_out, "cds"), os.path.join(s1_out, "gff"), os.path.join(s1_out, "bed")]

        if not step_done(s1_done, s1_required):
            print(f"[{get_now()}] [{gid}] Doing Step 1: Prepare Data...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "prepare_jcvi.py"),
                     "--gff", g_abs, "--fasta", f_abs, "--chr", str(chr_num), "--outdir", s1_out, "--protein"]):
                raise RuntimeError("Step 1 Failed")
            open(s1_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 1 Completed")

        # Step 1.5: Intra-Group OrthoFinder
        # CORE FIX: Removed the external .step1_5_done lock to allow the robust internal resume logic
        # of auto_of.py to handle the skip. This guarantees 'og_abs' is ALWAYS populated correctly.
        if g_args.auto_og:
            print(f"[{get_now()}] [{gid}] Doing Step 1.5: Intra-Group OrthoFinder (or checking cache)...")
            if PROJECT_ROOT not in sys.path:
                sys.path.insert(0, PROJECT_ROOT)
            from bin.auto_of import run_intra_group_orthofinder

            group_pep_files = glob.glob(os.path.join(s1_out, "**", "*.pep"), recursive=True)
            group_gff_files = glob.glob(os.path.join(s1_out, "**", "*.gff*"), recursive=True)

            if not group_pep_files or not group_gff_files:
                raise RuntimeError(f"Step 1.5 Failed: Cannot find .pep or .gff in {s1_out}")

            # Will intelligently resume or run from scratch, safely returning the required absolute path
            og_abs = run_intra_group_orthofinder(group_pep_files, g_dir, gid, g_args.threads)

        # Step 3: Tandem Duplication
        s3_out = os.path.join(g_dir, "02.tandem")
        m_pep = os.path.join(s3_out, f"{gid}_merged.pep")
        m_gff = os.path.join(s3_out, f"{gid}_merged.gff")
        tandem_file = os.path.join(s3_out, gid + ".tandem")
        s3_done = os.path.join(g_dir, ".step3_done")

        if not step_done(s3_done, [m_pep, m_gff, tandem_file], allow_empty_paths=[tandem_file]):
            print(f"[{get_now()}] [{gid}] Doing Step 3: Tandem Duplication Identification...")
            os.makedirs(s3_out, exist_ok=True)

            # Prevent FASTA/GFF format corruption by ensuring newline after each file
            with open(m_pep, 'wb') as w_pep, open(m_gff, 'wb') as w_gff:
                for pf in sorted(glob.glob(os.path.join(s1_out, "cds", "*.pep"))):
                    with open(pf, 'rb') as r:
                        shutil.copyfileobj(r, w_pep)
                    w_pep.write(b'\n')
                for gf in sorted(glob.glob(os.path.join(s1_out, "gff", "*.gff"))):
                    with open(gf, 'rb') as r:
                        shutil.copyfileobj(r, w_gff)
                    w_gff.write(b'\n')

            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "tandem_identify.py"),
                            "--gff", m_gff, "--pep", m_pep, "--outdir", s3_out, "--prefix", gid,
                            "--max_dist", str(g_args.tandem_dist),
                            "--threads", str(g_args.threads)]):
                raise RuntimeError("Step 3 Failed")
            open(s3_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 3 Completed")

        # Step 4: JCVI Synteny
        s4_out = os.path.join(g_dir, "03.jcvi")
        s4_done = os.path.join(g_dir, ".step4_done")

        if not step_done_with_glob(s4_done, [s4_out], [os.path.join(s4_out, "*.anchors")]):
            print(f"[{get_now()}] [{gid}] Doing Step 4: JCVI Synteny...")
            os.makedirs(s4_out, exist_ok=True)
            cds_src = os.path.join(s1_out, "cds")
            bed_src = os.path.join(s1_out, "bed")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "jcvi_anchors.py"),
                     "--cds_dir", cds_src, "--bed_dir", bed_src, "--jcvi_input", s4_out,
                     "--sh_script", os.path.join(BIN_DIR, "run_ortholog.sh")]):
                raise RuntimeError("Step 4 Failed")
            open(s4_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 4 Completed")

        # Step 5: Graph Clustering
        s5_out = os.path.join(g_dir, "04.cluster.tsv")
        s5_done = os.path.join(g_dir, ".step5_done")

        if not step_done(s5_done, [s5_out]):
            print(f"[{get_now()}] [{gid}] Doing Step 5: Graph Clustering...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "cluster.py"),
                     "--gff", m_gff, "--jcvi_dir", s4_out, "--tandem", tandem_file,
                     "--output", s5_out, "--min_chroms", str(g_args.min_c),
                     "--sub_list", g_args.sub_list, "--chr_num", str(chr_num),
                     "--max_gene_distance", str(g_args.cluster_dist)]):
                raise RuntimeError("Step 5 Failed")
            open(s5_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 5 Completed")

        # Step 6: Orthogonal Verification
        s6_pre = os.path.join(g_dir, "05.verified")
        s6_out = f"{s6_pre}.tsv"
        s6_done = os.path.join(g_dir, ".step6_done")

        if not step_done(s6_done, [s6_out]):
            print(f"[{get_now()}] [{gid}] Doing Step 6: Orthogonal Verification...")
            require_file(og_abs, "Orthogroups file")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "verify.py"),
                            "--cluster_file", s5_out, "--orthogroups", og_abs,
                            "--output_prefix", s6_pre, "--min_support_ratio", str(g_args.verify_ratio)]):
                raise RuntimeError("Step 6 Failed")
            open(s6_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 6 Completed")

        # Step 7: Sequence Homology Expansion
        s7_pre = os.path.join(g_dir, "06.expanded")
        s7_out = f"{s7_pre}_expanded.tsv"
        s7_done = os.path.join(g_dir, ".step7_done")

        if not step_done(s7_done, [s7_out]):
            print(f"[{get_now()}] [{gid}] Doing Step 7: Sequence Homology Expansion...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "blastout.py"),
                            "--allele_file", s6_out, "--fasta", m_pep,
                            "--out_prefix", s7_pre,
                            "--threads", str(g_args.threads)]):
                raise RuntimeError("Step 7 Failed")
            open(s7_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 7 Completed")

        # Step 8: Reference Calibration
        s8_out = os.path.join(g_dir, "07.FINAL_ALLELE.tsv")
        s8_done = os.path.join(g_dir, ".step8_done")

        if not step_done(s8_done, [s8_out]):
            print(f"[{get_now()}] [{gid}] Doing Step 8: Reference Calibration And Output...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "compare.py"),
                     "--allele_file", s7_out, "--ref_gff", ref_g_abs, "--ref_cds", ref_f_abs,
                     "--hap_cds", m_pep, "--output", s8_out, "--sub_list", g_args.sub_list,
                     "--chr_num", str(chr_num),
                     "--threads", str(g_args.threads)]):
                raise RuntimeError("Step 8 Failed")
            open(s8_done, 'w').close()
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 8 Completed")

        print(f"[{get_now()}] [SUCCESS] {gid} All steps completed strictly!")
        return True
    except Exception as e:
        print(f"[{get_now()}] [FAILED] {gid} Process Interruption: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="GraphAllele: Graph-Constrained Standardized Allele Matrix.")
    parser.add_argument("-g", "--gff", required=True)
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-ref_g", "--ref_gff", required=True)
    parser.add_argument("-ref_f", "--ref_cds", required=True)
    parser.add_argument("-og", "--orthogroups")
    parser.add_argument("--auto_og", action="store_true")
    parser.add_argument("-s", "--start", type=int, default=1)
    parser.add_argument("-e", "--end", type=int, default=10)
    parser.add_argument("-t", "--threads", type=int, default=10)
    parser.add_argument("-o", "--outdir", default="standardized_results")
    parser.add_argument("--min_c", type=int, default=3)
    parser.add_argument("--sub_list", default="A,B,C,D,E,F,G,H,I,J,K,L,M,N",
                        help="Preset haplotype suffix list for fixing column order")
    parser.add_argument("--tandem_dist", type=int, default=5,
                        help="Maximum gene index distance for tandem duplicates (MCScanX default: 5)")
    parser.add_argument("--cluster_dist", type=int, default=30,
                        help="Maximum gene index distance for synteny clustering (default: 30)")
    parser.add_argument("--verify_ratio", type=float, default=0.35,
                        help="Minimum Orthogroup support ratio for cluster verification (default: 0.35)")

    args = parser.parse_args()

    for path, label in [
        (args.gff, "input GFF"),
        (args.fasta, "input FASTA"),
        (args.ref_gff, "reference GFF"),
        (args.ref_cds, "reference CDS"),
    ]:
        require_file(path, label)

    if not args.auto_og:
        if not args.orthogroups:
            parser.error("Either --auto_og or -og/--orthogroups must be provided.")
        require_file(args.orthogroups, "orthogroups file")

    os.makedirs(args.outdir, exist_ok=True)

    task_list = [(i, args) for i in range(args.start, args.end + 1)]

    print(f"\n[{get_now()}] [INFO] Starting pipeline sequentially to protect HPC NFS I/O...")
    for task in task_list:
        success = run_group_pipeline(task)
        if not success:
            print(f"[{get_now()}] [FATAL ERROR] Pipeline aborted due to failure in Group_Chr{task[0]:02d}.")
            print(f"[{get_now()}] [ACTION] Please check logs, fix the issue, and rerun to auto-resume.")
            sys.exit(1)

    expected_files = []
    for i in range(args.start, args.end + 1):
        target_file = os.path.join(args.outdir, f"Group_Chr{i:02d}", "07.FINAL_ALLELE.tsv")
        if os.path.exists(target_file):
            expected_files.append(target_file)

    if expected_files:
        param_tag = f"minc{args.min_c}_dist{args.cluster_dist}_ratio{args.verify_ratio:g}"
        generate_clean_clusters_auto(
            work_dir=args.outdir,
            sub_list_str=args.sub_list,
            out_dense=f"PolyAlleler_Global_Matrix_Cleaned_{param_tag}.tsv",
            out_kmer=f"my_clusters_{param_tag}.tsv"
        )
    else:
        print(f"\n[{get_now()}] [WARNING] Failed to find any finalized TSV files for merging.")

if __name__ == "__main__":
    main()
