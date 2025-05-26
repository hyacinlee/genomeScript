#!/usr/bin/env python3
import re
import sys
from collections import defaultdict

def process_gff(gff_file, prefix, output):
    # 存储每个染色体的基因计数
    chr_gene_count = defaultdict(int)
    # 存储旧ID到新ID的映射
    id_mapping = {}
    
    # 第一遍扫描：统计每个染色体的基因数量并记录旧ID
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]
            
            if feature_type == "mRNA":
                # 提取染色体ID
                chr_id = parts[0]
                # 提取旧mRNA ID
                if 'ID=' in attributes:
                    old_id = re.search(r'ID=([^;]+)', attributes).group(1)
                    chr_gene_count[chr_id] += 1
    
    # 确定需要补零的长度
    max_count = max(chr_gene_count.values()) if chr_gene_count else 0
    zero_pad = len(str(max_count))
    
    # 重置计数器
    chr_gene_count = defaultdict(int)
    
    # 第二遍扫描：处理文件并输出
    with open(gff_file, 'r') as f_in, open(f"{output}.gff", 'w') as f_out, open(f"{output}.idchange", 'w') as f_out2:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                f_out.write(line)
                continue
            
            feature_type = parts[2]
            attributes = parts[8]
            
            if feature_type in ["mRNA", "CDS"]:
                # 提取染色体ID
                chr_id = parts[0]
                
                if feature_type == "mRNA":
                    if 'ID=' in attributes:
                        old_id = re.search(r'ID=([^;]+)', attributes).group(1)
                        chr_gene_count[chr_id] += 1
                        i = chr_gene_count[chr_id]
                        # 生成新ID
                        new_id = f"{prefix}0{chr_id}0G0{str(i).zfill(zero_pad)}"
                        f_out2.write(f"{new_id}\t{old_id}\n")
                        id_mapping[old_id] = new_id
                        # 替换ID
                        attributes = re.sub(r'ID=([^;]+)', f'ID={new_id}', attributes)
                elif feature_type == "CDS":
                    if 'Parent=' in attributes:
                        old_parent = re.search(r'Parent=([^;]+)', attributes).group(1)
                        if old_parent in id_mapping:
                            # 替换Parent
                            attributes = re.sub(r'Parent=([^;]+)', f'Parent={id_mapping[old_parent]}', attributes)
                
                # 更新属性字段
                parts[8] = attributes
                line = '\t'.join(parts) + '\n'
            
            f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.gff> <prefix> < output.gff>")
        sys.exit(1)
    
    gff_file = sys.argv[1]
    prefix = sys.argv[2]
    output = sys.argv[3]
    process_gff(gff_file, prefix, output)
    print(f"Processed map file saved as {output}.idchange")
    print(f"Processed GFF file saved as {output}.gff")
