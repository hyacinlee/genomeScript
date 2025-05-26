#!/usr/bin/env python3

#!/usr/bin/env python3
import re
from collections import defaultdict
import sys

def process_fasta(input_file, output_fasta, output_coords):
    # 存储所有contig信息：key为原始序列名，value为(contig_seq, start, end)列表
    all_contigs = defaultdict(list)
    
    # 读取并解析FASTA文件
    with open(input_file, 'r') as f:
        current_seq = None
        seq_data = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq is not None:
                    process_sequence(current_seq, ''.join(seq_data), all_contigs)
                current_seq = line[1:].split()[0]  # 只取头部第一个单词作为序列名
                seq_data = []
            else:
                seq_data.append(line)
        
        # 处理最后一个序列
        if current_seq is not None:
            process_sequence(current_seq, ''.join(seq_data), all_contigs)
    
    # 收集所有contig并排序
    sorted_contigs = []
    for seq_name, contigs in all_contigs.items():
        for contig in contigs:
            # 现在contig的结构是 (start, end, sequence)
            sorted_contigs.append((seq_name, contig[0], contig[1], contig[2]))
    
    # 按长度降序排序 (现在x[3]是序列字符串)
    sorted_contigs.sort(key=lambda x: len(x[3]), reverse=True)
    
    # 写入新的FASTA文件和坐标文件
    with open(output_fasta, 'w') as fasta_out, open(output_coords, 'w') as coords_out:
        coords_out.write("#ContigID\tCstart\tLength\tOriginalSeq\tStart\tEnd\n")
        
        for i, (seq_name, start, end, contig) in enumerate(sorted_contigs, 1):
            contig_id = f"contig{i}"
            length = len(contig)
            
            # 写入FASTA
            fasta_out.write(f">{contig_id}\n")
            for j in range(0, len(contig), 80):
                fasta_out.write(contig[j:j+80] + '\n')
            
            # 写入坐标
            coords_out.write(f"{contig_id}\t1\t{length}\t{seq_name}\t{start}\t{end}\n")

def process_sequence(seq_name, sequence, all_contigs):
    # 找到所有非N的连续区域
    non_n_regions = [(m.start(), m.end()-1) for m in re.finditer(r'[^Nn]+', sequence)]
    
    for start, end in non_n_regions:
        contig_seq = sequence[start:end+1]
        all_contigs[seq_name].append((start+1, end+1, contig_seq))  # 转换为1-based坐标

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.fasta> <output.fasta> <output.coords>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_fasta = sys.argv[2]
    output_coords = sys.argv[3]
    
    process_fasta(input_file, output_fasta, output_coords)
    print(f"Processing complete. Contigs written to {output_fasta}, coordinates to {output_coords}")
