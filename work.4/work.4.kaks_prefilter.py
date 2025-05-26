#!/usr/bin/env python3
import sys
import re

# 假设文件名为 'genes.txt'
filename = sys.argv[1]
out   = sys.argv[2]
# 定义一个函数来提取基因名中的 'n' 值
def extract_n(gene_name):
    # 正则表达式匹配 xxx0xx{n}0Gxxxx 中的 n
    match = re.search(r'(\d)0G', gene_name)
    if match:
        return match.group(1)  # 返回 n 的值
    return None

# 打开文件并逐行读取
with open(filename, 'r') as f, open(f"{out}.A_B.lst","w") as fab, open(f"{out}.A_C.lst","w") as fac,open(f"{out}.B_C.lst","w") as fbc:
    for line in f:
        # 以空格或制表符为分隔符，获取各列
        columns = line.strip().split()
        # 获取每行中从第二列开始的基因名
        genes = columns
        # 提取每个基因名中的 'n' 并检查它们是否相同
        n_values = [extract_n(gene) for gene in genes]
        
        # 如果所有 n 值相同
        if len(set(n_values)) == 1 and n_values[0] is not None:
            fab.write(f"{genes[0]}\t{genes[1]}\t{genes[2]}\n")
            fac.write(f"{genes[0]}\t{genes[1]}\t{genes[3]}\n")
            fbc.write(f"{genes[0]}\t{genes[2]}\t{genes[3]}\n")
            continue
            print(f"行 '{line.strip()}' 中的基因 n 值相同: {n_values[0]}")
        else:
            print(f"行 '{line.strip()}' 中的基因 n 值不同或无效")

