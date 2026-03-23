# -*- coding: utf-8 -*-
"""
PolyAlleler V6.2 - Standardized Allele Matrix Version (Graph-Constrained)
Function: Supports preset subgenome structures, breakpoint resume, and dynamic graph constraints.
Fixed: Sequential processing to prevent HPC NFS I/O deadlock, early import validation, and strict error handling.
"""

import os
import sys
import subprocess
import argparse
import shutil
import datetime
import glob
import pandas as pd  # 【修复2】提前暴露环境问题：把 pandas 移到文件头部，没有装立刻报错，绝不浪费30小时机时

WORKFLOW_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(WORKFLOW_DIR)
BIN_DIR = os.path.join(PROJECT_ROOT, "bin")

def get_now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def run_cmd(cmd, cwd=None):
    try:
        subprocess.run(cmd, check=True, cwd=cwd)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[{get_now()}] [ERROR] Command failed: {' '.join(cmd)}")
        return False

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
        # --- Step 1: Prepare Data ---
        s1_out = os.path.join(g_dir, "01.prepare")
        s1_check = os.path.join(s1_out, "cds")
        if not os.path.exists(s1_check):
            print(f"[{get_now()}] [{gid}] Doing Step 1: Prepare Data...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "prepare_jcvi.py"), 
                     "--gff", g_abs, "--fasta", f_abs, "--chr", str(chr_num), "--outdir", s1_out, "--protein"]):
                raise RuntimeError("Step 1 Failed")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 1 Completed, automatically skip")

        # ==================== Step 1.5 Intra-Group OrthoFinder ====================
        if g_args.auto_og:
            print(f"[{get_now()}] [{gid}] Doing Step 1.5: Intra-Group OrthoFinder...")
            sys.path.append(PROJECT_ROOT)
            from bin.auto_of import run_intra_group_orthofinder
            
            group_pep_files = glob.glob(os.path.join(s1_out, "**", "*.pep"), recursive=True)
            group_gff_files = glob.glob(os.path.join(s1_out, "**", "*.gff*"), recursive=True)
            
            if not group_pep_files or not group_gff_files:
                raise RuntimeError(f"Step 1.5 Failed: Cannot find .pep or .gff in {s1_out}")
                
            # 【核心护城河】：此时 g_args.threads 会把PBS申请的全部线程（如20核）毫无保留地砸给当前这唯一的 Group，极速狂飙！
            og_abs = run_intra_group_orthofinder(group_pep_files, g_dir, gid, g_args.threads)
        # ===============================================================================

        # --- Step 3: Tandem ---
        s3_out = os.path.join(g_dir, "02.tandem")
        m_pep = os.path.join(s3_out, f"{gid}_merged.pep")
        m_gff = os.path.join(s3_out, f"{gid}_merged.gff")
        tandem_file = os.path.join(s3_out, gid + ".tandem")
        
        if not os.path.exists(tandem_file):
            print(f"[{get_now()}] [{gid}] Doing Step 3: Tandem Duplication Identification...")
            os.makedirs(s3_out, exist_ok=True)
            with open(m_pep, 'wb') as w_pep, open(m_gff, 'wb') as w_gff:
                for pf in sorted(glob.glob(os.path.join(s1_out, "cds", "*.pep"))):
                    with open(pf, 'rb') as r: shutil.copyfileobj(r, w_pep)
                for gf in sorted(glob.glob(os.path.join(s1_out, "gff", "*.gff"))):
                    with open(gf, 'rb') as r: shutil.copyfileobj(r, w_gff)
            
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "tandem_identify.py"), 
                            "--gff", m_gff, "--pep", m_pep, "--outdir", s3_out, "--prefix", gid,
                            "--max_dist", str(g_args.tandem_dist)]):
                raise RuntimeError("Step 3 Failed")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 3 Completed, automatically skip")

        # --- Step 4: JCVI ---
        s4_out = os.path.join(g_dir, "03.jcvi")
        jcvi_done_mark = os.path.join(s4_out, ".jcvi_done")
        if not os.path.exists(jcvi_done_mark):
            print(f"[{get_now()}] [{gid}] Doing Step 4: JCVI Synteny...")
            os.makedirs(s4_out, exist_ok=True)
            cds_src = os.path.join(s1_out, "cds")
            bed_src = os.path.join(s1_out, "bed")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "jcvi_anchors.py"), 
                     "--cds_dir", cds_src, "--bed_dir", bed_src, "--jcvi_input", s4_out, 
                     "--sh_script", os.path.join(BIN_DIR, "run_ortholog.sh")]):
                raise RuntimeError("Step 4 Failed")
            with open(jcvi_done_mark, 'w') as f: f.write("done")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 4 Completed, automatically skip")

        # --- Step 5: Graph Clustering ---
        s5_out = os.path.join(g_dir, "04.cluster.tsv")
        if not os.path.exists(s5_out):
            print(f"[{get_now()}] [{gid}] Doing Step 5: Graph Clustering...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "cluster.py"),
                     "--gff", m_gff, "--jcvi_dir", s4_out, "--tandem", tandem_file,
                     "--output", s5_out, "--min_chroms", str(g_args.min_c), 
                     "--sub_list", g_args.sub_list, "--chr_num", str(chr_num),
                     "--max_gene_distance", str(g_args.cluster_dist)]):
                raise RuntimeError("Step 5 Failed")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 5 Completed, automatically skip")

        # --- Step 6: Verify ---
        s6_pre = os.path.join(g_dir, "05.verified")
        s6_out = f"{s6_pre}.tsv"
        if not os.path.exists(s6_out):
            print(f"[{get_now()}] [{gid}] Doing Step 6: Verify...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "verify.py"), "--cluster_file", s5_out, "--orthogroups", og_abs, "--output_prefix", s6_pre]):
                raise RuntimeError("Step 6 Failed")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 6 Completed, automatically skip")

        # --- Step 7: Expanded ---
        s7_pre = os.path.join(g_dir, "06.expanded")
        s7_out = f"{s7_pre}_expanded.tsv"
        if not os.path.exists(s7_out):
            print(f"[{get_now()}] [{gid}] Doing Step 7: Sequence Homology Expansion...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "blastout.py"), "--allele_file", s6_out, "--fasta", m_pep, "--out_prefix", s7_pre]):
                raise RuntimeError("Step 7 Failed")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 7 Completed, automatically skip")

        # --- Step 8: Reference Calibration ---
        s8_out = os.path.join(g_dir, "07.FINAL_ALLELE.tsv")
        if not os.path.exists(s8_out):
            print(f"[{get_now()}] [{gid}] Doing Step 8: Reference Calibration And Output...")
            if not run_cmd([sys.executable, os.path.join(BIN_DIR, "compare.py"),
                     "--allele_file", s7_out, "--ref_gff", ref_g_abs, "--ref_cds", ref_f_abs,
                     "--hap_cds", m_pep, "--output", s8_out, "--sub_list", g_args.sub_list, "--chr_num", str(chr_num)]):
                raise RuntimeError("Step 8 Failed")
        else:
            print(f"[{get_now()}] [SKIP] {gid} Step 8 Completed, automatically skip")

        print(f"[{get_now()}] [SUCCESS] {gid} All steps completed!")
        return True
    except Exception as e:
        print(f"[{get_now()}] [FAILED] {gid} Process Interruption: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="PolyAlleler V6.2: Graph-Constrained Standardized Allele Matrix. Developed by engineers Zhangqing-Lab, Yi-Chen, and Gengrui Zhu. We thank them for their outstanding contributions to polyploid research.")
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

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    task_list = [(i, args) for i in range(args.start, args.end + 1)]

    # 【修复1】彻底废弃 Pool 并发，改为稳如泰山的串行（Sequential）运行！
    # 彻底杜绝超算网络硬盘（NFS）的 I/O 死锁踩踏惨案。
    print(f"\n[{get_now()}] [INFO] Starting pipeline sequentially to protect HPC NFS I/O...")
    for task in task_list:
        success = run_group_pipeline(task)
        # 【修复3】严格报错机制：一旦某一个群彻底失败，直接终止主程序，绝不带着残缺数据往后跑！
        if not success:
            print(f"[{get_now()}] [FATAL ERROR] Pipeline aborted due to failure in Group_Chr{task[0]:02d}.")
            print(f"[{get_now()}] [ACTION] Please check the error logs above, fix the issue, and rerun. The script will automatically resume from this group.")
            sys.exit(1)

    expected_files = []
    for i in range(args.start, args.end + 1):
        target_file = os.path.join(args.outdir, f"Group_Chr{i:02d}", "07.FINAL_ALLELE.tsv")
        if os.path.exists(target_file):
            expected_files.append(target_file)

    if len(expected_files) > 1:
        print(f"\n[{get_now()}] [INFO] Testing having more than two Groups, starting Merging...")
        df_list = [pd.read_csv(f, sep='\t') for f in expected_files]
        global_df = pd.concat(df_list, ignore_index=True)
        
        # 重新编号 ClusterID，打造全局唯一的家族 ID
        global_df['ClusterID'] = [f"Global_Cluster_{i:06d}" for i in range(len(global_df))]
        
        global_out = os.path.join(args.outdir, "PolyAlleler_Global_Matrix.tsv")
        global_df.to_csv(global_out, sep='\t', index=False)
        
        print(f"[{get_now()}] [SUCCESS] Final tsv: {global_out}")
        print(f"[{get_now()}] [SUCCESS] Allele groups: {len(global_df)} 个。")
    elif len(expected_files) == 1:
        print(f"\n[{get_now()}] [INFO] Just Single Group, skipping merging.")
    else:
        print(f"\n[{get_now()}] [WARNING] Failed to find any 07.FINAL_ALLELE.tsv files.")

if __name__ == "__main__":
    main()
