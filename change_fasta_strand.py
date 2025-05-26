#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def process_fasta_with_table(fasta_file, table_file, output_fasta, output_table):
    """
    根据表格重命名 FASTA 并调整链方向，输出新 FASTA 和表格。
    """
    # 读取表格
    df = pd.read_csv(table_file, sep='\t', header=None, names=['new_id', 'old_id', 'strand'])
    
    # 读取 FASTA 并建立旧 ID 到序列的映射
    fasta_records = {rec.id: rec for rec in SeqIO.parse(fasta_file, 'fasta')}
    
    # 处理并输出
    with open(output_fasta, 'w') as out_fasta, open(output_table, 'w') as out_table:
        out_table.write("#new_id\told_id\tStart\tEnd\tstrand\n")  # 表头
        
        for _, row in df.iterrows():
            old_id = row['old_id']
            new_id = row['new_id']
            strand = row['strand']
            
            if old_id not in fasta_records:
                print(f"Warning: {old_id} not found in FASTA, skipping.", file=sys.stderr)
                continue
            
            # 获取序列并调整方向
            record = fasta_records[old_id]
            seq = record.seq
            if strand == '-':
                seq = seq.reverse_complement()
            
            # 写入新 FASTA
            new_record = SeqRecord(
                seq,
                id=new_id,
                description=f"old_id={old_id} strand={strand}"
            )
            SeqIO.write(new_record, out_fasta, 'fasta')
            
            # 写入表格
            out_table.write(f"{new_id}\t{old_id}\t1\t{len(seq)}\t{strand}\n")

process_fasta_with_table(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
