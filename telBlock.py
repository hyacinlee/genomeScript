#!/usr/bin/env python 
import sys 
import re 

def main():

    if len(sys.argv) < 3 :
        print "# <Usage>  xxx.fasta xxx.fasta.telomere  "
        exit(1)

    fasta=readFasta(sys.argv[1])
    for line in open(sys.argv[2]):
        if line.startswith("#"):
            continue
        info=line.strip().split()
        sb="%s\t%s\t%s" % (info[0],1,10000)
        eb="%s\t%s\t%s" % (info[0],fasta[info[0]]-10000,fasta[info[0]])
        if info[6] == "Both":
            print sb 
            print eb
        elif info[6] == "Single-start":
            print sb 
        elif info[6] == "Single-end":
            print eb 
        else:
            continue



def readFasta(fasta_file):
    fasta    = {}
    fasta_id = ''
    for line in open(fasta_file):
        if line.startswith(">"):
            fasta_id = line.strip().replace(">","")
            fasta[fasta_id] = 0
        else:
            fasta[fasta_id] += len(line.strip())
    return fasta 


main()
