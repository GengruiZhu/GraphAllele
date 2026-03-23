# -*- coding: utf-8 -*-
"""
GraphAllele 稀疏矩阵降维清洗器
功能：
1. 将 140 列的稀疏多倍体矩阵，压缩为标准的 A-N 14列亚基因组密集矩阵。
2. 自动生成下游 EM 定量引擎所需的 clusters_table.tsv。
"""
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("用法: python clean_allele_matrix.py <PolyAlleler_Global_Matrix.tsv>")
        sys.exit(1)

    input_file = sys.argv[1]
    
    # 自动定义输出文件名
    base_name = os.path.splitext(input_file)[0]
    dense_out = f"{base_name}_Cleaned.tsv"       # 给你看的完美表格
    cluster_out = "my_clusters.tsv"              # 给定量引擎吃的配置文件

    print(f"[INFO] 正在读取原始稀疏矩阵: {input_file}")

    with open(input_file, 'r') as fin, \
         open(dense_out, 'w') as fout_dense, \
         open(cluster_out, 'w') as fout_cluster:

        # 1. 解析表头，建立“列索引 -> 亚基因组(A-N)”的映射雷达
        header = fin.readline().strip('\n').split('\t')
        sg_map = {}
        for i in range(3, len(header)):
            col = header[i].strip()
            # 匹配形如 Chr01A, Chr10N 的列名，提取最后一个字母
            if len(col) >= 6 and col.startswith("Chr"):
                sg = col[-1]
                if 'A' <= sg <= 'N':
                    sg_map[i] = sg

        # 强制按 A 到 N 的顺序生成 14 个槽位
        subgenomes = [chr(i) for i in range(ord('A'), ord('N')+1)]

        # 写入 Cleaned 矩阵表头
        dense_header = ['ClusterID', 'Ref_Gene', 'Ref_Locus'] + [f'Allele_{sg}' for sg in subgenomes]
        fout_dense.write('\t'.join(dense_header) + '\n')

        line_count = 0
        for line in fin:
            if not line.strip(): continue
            parts = line.strip('\n').split('\t')
            if len(parts) < 3: continue

            cid = parts[0].strip()
            ref_g = parts[1].strip()
            ref_l = parts[2].strip()

            # 初始化该 Cluster 的 14 个亚基因组槽位 (默认填横杠)
            allele_dict = {sg: '-' for sg in subgenomes}
            valid_alleles = [] # 收集非空的有效基因 ID

            # 遍历该行的所有列（从第4列开始）
            for i in range(3, len(parts)):
                val = parts[i].strip()
                if val and val != '-':
                    if i in sg_map:
                        sg = sg_map[i]
                        # 填入对应的亚基因组槽位 (如果同一个槽位碰到多个基因，用逗号连起来，防止串联重复被覆盖)
                        if allele_dict[sg] == '-':
                            allele_dict[sg] = val
                        else:
                            allele_dict[sg] += f",{val}"
                        
                        valid_alleles.append(val)

            # --- 动作 1：写入人类可读的压缩矩阵 ---
            row = [cid, ref_g, ref_l] + [allele_dict[sg] for sg in subgenomes]
            fout_dense.write('\t'.join(row) + '\n')

            # --- 动作 2：生成给引擎吃的输入文件 (只保留 ClusterID 和用逗号分隔的基因列表) ---
            if valid_alleles:
                fout_cluster.write(f"{cid}\t{','.join(valid_alleles)}\n")
            
            line_count += 1

    print(f"[SUCCESS] 矩阵降维完成！共处理 {line_count} 个 Clusters。")
    print(f"  -> 人类可读报表已保存至: {dense_out}")
    print(f"  -> 定量引擎配置已保存至: {cluster_out}")

if __name__ == "__main__":
    main()
