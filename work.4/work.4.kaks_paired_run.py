#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from multiprocessing import Pool
import shutil

def process_row(row, pfasta_dict,cfasta_dict,intab_tmp_dir):
    name, gene1, gene2 = row
    output_file = os.path.join(intab_tmp_dir, f"{name}")
    poutput_file = os.path.join(intab_tmp_dir, f"{name}.pep.fasta")
    coutput_file = os.path.join(intab_tmp_dir, f"{name}.cds.fasta")

    # 获取序列并写入临时文件
    with open(poutput_file, "w") as fp,open(coutput_file, "w") as fc:
        pfa1="".join(pfasta_dict[gene1])
        cfa1="".join(cfasta_dict[gene1])
        pfa2="".join(pfasta_dict[gene2])
        cfa2="".join(cfasta_dict[gene2])

        fp.write(f">{gene1}\n{pfa1}\n>{gene2}\n{pfa2}")
        fc.write(f">{gene1}\n{cfa1}\n>{gene2}\n{cfa2}")


    os.system(f"mafft  {poutput_file} 1>{poutput_file}.mafft 2>{poutput_file}.e")
    os.system(f"Epal2nal.pl {poutput_file}.mafft {coutput_file} -nogap -output fasta > {coutput_file}.conda.aln")

    axt_file=os.path.join(intab_tmp_dir, f"{name}.axt")
    create_axt_file(f"{coutput_file}.conda.aln",axt_file,name)

    os.system(f"KaKs_Calculator -i  {axt_file} -o {axt_file}.kaks -m GMYN")
    os.system(f"rm -rf {output_file}.cds.* {output_file}.pep.*")

    return f'{axt_file}.kaks'

def main(pfasta_file,cfasta_file,intable_file, threads=32):
    # 读取FASTA文件
    pfasta_dict = read_fasta(pfasta_file)
    cfasta_dict = read_fasta(cfasta_file)
    # 读取表格文件
    infs=[]
    for l in open(intable_file,"r"):
        ls = l.strip().split("\t")
        infs.append(ls)
    
    # 创建临时目录
    intab_tmp_dir = f"{os.path.splitext(intable_file)[0]}.tmp"
    os.makedirs(intab_tmp_dir, exist_ok=True)
    
    # 准备参数列表
    tasks = [(inf, pfasta_dict,cfasta_dict,intab_tmp_dir) for inf in infs]
    
    # 使用多进程并行处理
    with Pool(processes=threads) as pool:
        results = pool.starmap(process_row, tasks)
    
    # 打印结果统计
    successful = sum(1 for r in results if r is not None)
    print(f"\nProcessing completed. Success: {successful}/{len(infs)}")
    
    head_flag=0
    with open(f"{intable_file}.kaks.result","w") as fo:
        for f in results:
            with open(f,"r") as fi:
                lines = [line.strip() for line in fi]
                if head_flag ==0:
                    fo.write(f"{lines[0]}\n")
                    head_flag = 1 
                if len(lines) > 1:
                    fo.write(f"{lines[1]}\n")
   

    # 返回临时目录路径，方便后续操作
    return intab_tmp_dir

def create_axt_file(fasta_file, output_file,name):
    # 读取FASTA文件
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) != 2:
        raise ValueError("FASTA文件应包含且仅包含2条序列")
    
    # 获取序列并去除可能的空白字符
    seq1 = str(records[0].seq).strip().upper()
    seq2 = str(records[1].seq).strip().upper()
    
    # 验证长度
    if len(seq1) != len(seq2):
        raise ValueError(f"序列长度不一致: {len(seq1)} vs {len(seq2)}")
    
    # 验证是否为有效DNA序列
    valid_bases = {'A', 'T', 'C', 'G'}
    if any(base not in valid_bases for base in seq1+seq2):
        raise ValueError("序列包含非标准碱基")
    
    # 写入AXT文件
    with open(output_file, 'w') as f:
        f.write(f"{name}\n")
        f.write(f"{seq1}\n")
        f.write(f"{seq2}\n")



def read_fasta(fasta_file,order=None):
    '''
        reading fasta file and return a dict of name and sequence .
    '''
    #print "# Reading fasta dict from %s" %(fasta_file)
    fasta   = {}
    fasta_id = ''
    fo=[]
    for line in open(fasta_file):
        if line.startswith(">"):
            fasta_id = line.strip().replace(">","").split()[0]
            fasta[fasta_id] = []
            fo.append(fasta_id)
        else:
            fasta[fasta_id].append(line.strip())
    #print "# Finish reading %s seqs in %s" %(len(fasta.keys()),fasta_file)
    #print fo
    if not order:
        return fasta 
    else:
        return fasta,fo



if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py < pep fasta file > < cds fasta file > < intable file > [threads]")
        sys.exit(1)
    
    pfasta_file = sys.argv[1]
    cfasta_file = sys.argv[2]
    intable_file = sys.argv[3]
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 4

    print(sys.argv)
    
    if not os.path.exists(pfasta_file):
        print(f"Error: FASTA file {fasta_file} not found")
        sys.exit(1)
    
    if not os.path.exists(intable_file):
        print(f"Error: Table file {intable_file} not found")
        sys.exit(1)
    
    output_dir = main(pfasta_file, cfasta_file, intable_file, threads)
    print(f"Alignment files saved in: {output_dir}")
