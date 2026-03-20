#!/bin/bash
# 工作目录（含 .cds 和 .bed 文件）
INPUT_DIR="./"
CSCORE=0.99

# 自动提取物种名（前提：同名 .cds 和 .bed 都存在）
species=()
for cds_file in "$INPUT_DIR"/*.cds; do
    base=$(basename "$cds_file" .cds)
    bed_file="$INPUT_DIR/$base.bed"
    if [[ -f "$bed_file" ]]; then
        species+=("$base")
    fi
done

# 检查是否有足够物种
if [[ ${#species[@]} -lt 2 ]]; then
    echo "[ERROR] Less than 2 valid species found!"
    exit 1
fi

# 两两比对
for (( i=0; i<${#species[@]}; i++ )); do
    for (( j=i+1; j<${#species[@]}; j++ )); do
        a="${species[$i]}"
        b="${species[$j]}"
        echo "[INFO] Running JCVI ortholog: $a vs $b"
        # 【核心修改点】：在命令最末尾加上 --no_plot，强制禁止画图，避免 zlib 崩溃
        python -m jcvi.compara.catalog ortholog --no_strip_names "$a" "$b" --cscore=$CSCORE --no_dotplot
    done
done

echo "[DONE] JCVI ortholog analysis completed!"
