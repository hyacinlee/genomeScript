#!/usr/bin/env python3
import argparse
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def extract_number(chr_name):
    """从染色体名中提取数字部分（如 'chr1' -> 1，'scaffold_2' -> 2）"""
    match = re.search(r'\d+', chr_name)
    return int(match.group()) if match else 0  # 若无数字则默认为0

def process_coords_and_fasta(coords_file, b_fasta, output_fasta, mapping_file, chr_prefix):
    """
    功能：
    1. 仅输出有比对的 B 染色体。
    2. 重命名为 chr_prefix + i（i为A染色体数字部分）。
    3. 按数字升序排序。
    """
    best_mappings = defaultdict(lambda: {
        'ref_chr': None,
        'strand': '+',
        'len': 0,
        'cov_q': 0.0,
        'num': 0
    })
    used_b_chrs = set()  # 记录已比对的B染色体

    # 第一步：解析 coords 文件
    with open(coords_file, 'r') as f:
        for line in f:
            if line.startswith('[') or not line.strip():
                continue
            parts = re.split(r'\s+\|\s+|\s+', line.strip())
            s1, e1, s2, e2 = map(int, parts[:4])
            len1, len2 = map(int, parts[4:6])
            cov_q = float(parts[8])  # COV Q
            tag_a, tag_b = parts[-2], parts[-1]
            strand = '+' if s2 < e2 else '-'
            num = extract_number(tag_a)
            used_b_chrs.add(tag_b)  # 标记为已比对

            # 更新最佳匹配
            current_best = best_mappings[tag_b]
            if len2 > current_best['len'] or (len2 == current_best['len'] and cov_q > current_best['cov_q']):
                best_mappings[tag_b] = {
                    'ref_chr': tag_a,
                    'strand': strand,
                    'len': len2,
                    'cov_q': cov_q,
                    'num': num
                }

    # 第二步：输出前检查B染色体是否在比对中
    all_b_chrs = {rec.id for rec in SeqIO.parse(b_fasta, 'fasta')}
    unused_b_chrs = all_b_chrs - used_b_chrs
    if unused_b_chrs:
        print(f"Warning: {len(unused_b_chrs)} B chromosomes have no alignment and will be skipped:",
              ", ".join(sorted(unused_b_chrs)), file=sys.stderr)

    # 第三步：按数字升序输出有比对的染色体
    with open(output_fasta, 'w') as out_fasta, open(mapping_file, 'w') as out_map:
        out_map.write("#NewID\tOldID\tStart\tEnd\tStrand\tOriginalRef\n")
        
        # 按数字升序处理
        sorted_b_chrs = sorted(
            [b for b in best_mappings.keys() if best_mappings[b]['ref_chr']],  # 仅保留有ref_chr的
            key=lambda x: best_mappings[x]['num']
        )
        
        for b_chr in sorted_b_chrs:
            mapping = best_mappings[b_chr]
            seq_record = next(rec for rec in SeqIO.parse(b_fasta, 'fasta') if rec.id == b_chr)
            seq = seq_record.seq
            if mapping['strand'] == '-':
                seq = seq.reverse_complement()

            new_id = f"{chr_prefix}{mapping['num']}"
            new_record = SeqRecord(
                seq,
                id=new_id,
                description=f"original_id={b_chr} strand={mapping['strand']} ref={mapping['ref_chr']}"
            )
            SeqIO.write(new_record, out_fasta, 'fasta')
            out_map.write(f"{new_id}\t{b_chr}\t1\t{len(seq)}\t{mapping['strand']}\t{mapping['ref_chr']}\n")

def main():
    parser = argparse.ArgumentParser(description='Adjust B genome with filtering unaligned chromosomes.')
    parser.add_argument('--coords', required=True, help='MUMmer coords file')
    parser.add_argument('--b_fasta', required=True, help='B genome FASTA file')
    parser.add_argument('--output', required=True, help='Output adjusted FASTA file')
    parser.add_argument('--mapping', required=True, help='Output mapping file (TSV)')
    parser.add_argument('--chr_prefix', default='chr', help='Prefix for renamed chromosomes (e.g., "chr" or "scaffold")')
    args = parser.parse_args()

    process_coords_and_fasta(args.coords, args.b_fasta, args.output, args.mapping, args.chr_prefix)
    print(f"Adjusted B genome (aligned chromosomes only): {args.output}")
    print(f"Mapping table: {args.mapping}")

if __name__ == '__main__':
    main()
