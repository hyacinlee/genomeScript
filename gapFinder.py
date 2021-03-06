#!/usr/bin/env python 
import sys 
import re 

def main():

    if len(sys.argv) < 2 :
        print "# please give a fasta as input ! "
        exit(1)

    fasta=readFasta(sys.argv[1])
    for seq  in fasta.keys():
        for match in re.finditer(r'N+',"".join(fasta[seq])):
            print "%s\t%s\t%s\t%s" % (seq,match.start()-1, match.end()-2,match.end()-match.start()+1) 





def readFasta(fasta_file):
    #'''
    #reading fasta file and return a dict of name and sequence .
    #'''
    print "# Reading fasta dict from %s" %(fasta_file)
    fasta    = {}
    fasta_id = ''
    for line in open(fasta_file):
        if line.startswith(">"):
            fasta_id = line.strip().replace(">","")
            fasta[fasta_id] = []
        else:
            fasta[fasta_id].append(line.strip())
    #print "# Finish reading %s seqs in %s" %(len(fasta.keys()),fasta_file)
    return fasta 


main()
