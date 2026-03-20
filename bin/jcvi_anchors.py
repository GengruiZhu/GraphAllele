# -*- coding: utf-8 -*-
import os
import shutil
import subprocess
import argparse
import sys

def copy_files(src_dir, dst_dir, extensions):
    """
    将所需的序列和位置文件搬运到 JCVI 的工作目录下
    """
    os.makedirs(dst_dir, exist_ok=True)
    src_dir = os.path.abspath(src_dir)
    dst_dir = os.path.abspath(dst_dir)
    
    count = 0
    for fname in os.listdir(src_dir):
        if any(fname.endswith(ext) for ext in extensions):
            src = os.path.join(src_dir, fname)
            dst = os.path.join(dst_dir, fname)
            # 使用 symlink 替代 copy2，既节省空间又防止权限问题
            if os.path.exists(dst): os.remove(dst)
            os.symlink(src, dst)
            count += 1
    print(f"[INFO] Linked {count} files ({extensions}) to {dst_dir}")

def run_shell_script(script_path, work_dir):
    """
    协同核心：必须进入 work_dir 后再运行 Shell 脚本
    """
    script_abs_path = os.path.abspath(script_path)
    print(f"[INFO] 正在进入工作目录: {work_dir}")
    print(f"[INFO] 准备运行 JCVI 脚本: {script_abs_path}")
    
    # 关键：通过 cwd 参数确保 Shell 脚本在有文件的目录下运行
    result = subprocess.run(["bash", script_abs_path], cwd=work_dir)
    if result.returncode != 0:
        raise Exception(f"Shell 脚本执行失败，退出码: {result.returncode}")

def collect_anchors(input_dir, output_dir):
    """
    收割比对结果 (.anchors)
    """
    os.makedirs(output_dir, exist_ok=True)
    count = 0
    for root, _, files in os.walk(input_dir):
        for f in files:
            # 只要原始锚点，不要 lifted 后的
            if f.endswith(".anchors") and not f.endswith(".lifted.anchors"):
                src = os.path.join(root, f)
                dst = os.path.join(output_dir, f)
                shutil.copy2(src, dst)
                count += 1
    print(f"[INFO] 成功收割 {count} 个 .anchors 文件至 {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="JCVI Pipeline Orchestrator for PolyAlleler")
    parser.add_argument("--cds_dir", required=True, help="Directory containing .cds files")
    parser.add_argument("--bed_dir", required=True, help="Directory containing .bed files")
    parser.add_argument("--jcvi_input", default="jcvi_work", help="Working dir for JCVI process")
    parser.add_argument("--anchors_dir", default="anchors_output", help="Final anchors output dir")
    parser.add_argument("--sh_script", required=True, help="Path to run_ortholog.sh")
    args = parser.parse_args()

    # 1. 搬运原料
    # JCVI 要求 .cds 和 .bed 在同一目录下
    copy_files(args.cds_dir, args.jcvi_input, [".cds"])
    copy_files(args.bed_dir, args.jcvi_input, [".bed"])

    # 2. 运行比对 (进入工作区)
    try:
        run_shell_script(args.sh_script, args.jcvi_input)
    except Exception as e:
        print(f"[ERROR] JCVI 运行过程中发生崩溃: {e}")
        sys.exit(1)

    # 3. 收割结果
    collect_anchors(args.jcvi_input, args.anchors_dir)
    print("[SUCCESS] JCVI 环节任务圆满完成。")

if __name__ == "__main__":
    main()
