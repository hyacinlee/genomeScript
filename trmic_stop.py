#!/usr/bin/env python 
import sys 


def main():
    if not  len(sys.argv) == 4:
        print "%s\t < input cds > < input gff > < out.prefix >" %(sys.argv[0])
        exit(1)
    stops=["TAA","TGA","TAG"] #TAA' => 'U', 'TAG' => 'U', 'TGA
    fasta,fid = readFasta(sys.argv[1],1)
    for ids in fid:
        newseq=""
        #sid=0
        for tri in conda("".join(fasta[ids]),3):
            newseq+=tri
            if tri in stops:
               #sid=i
               break
        if not len(newseq) == len("".join(fasta[ids])): 
            print "%s\t%s\t%s\t%s" % (ids,len("".join(fasta[ids])),len(newseq),tri)


def conda(xxx,bp=3):
    if not len(xxx) % bp == 0:
        print "#Error: Con't be exact division by %s: %s" %(bp,xxx)
        exit(1)
    ll=[]
    for i in range(len(xxx)/bp):
        s= i*bp 
        e=(i+1)*bp
        ll.append(xxx[s:e])
    return ll

 
def readFasta(fasta_file,order=0):
    '''
        reading fasta file and return a dict of name and sequence .
    '''
    #print "# Reading fasta dict from %s" %(fasta_file)
    fo       = []
    fasta    = {}
    fasta_id = ''
    for line in open(fasta_file,"r"):
        if line.startswith(">"):
            fasta_id = line.strip().replace(">","").split()[0]
            fasta[fasta_id] = []
            fo.append(fasta_id)
        else:
            fasta[fasta_id].append(line.strip())
    #print "# Finish reading %s seqs in %s" %(len(fasta.keys()),fasta_file)
    #print fo
    if order == 0:
        return fasta 
    else:
        return fasta,fo




main()





