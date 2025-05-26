#!/usr/bin/env python 
import sys 
if len(sys.argv) < 2 :
    print "\tFunction: remove dup-busco in short seqs\n\tUsage: python rmdupBusco.py  < busco table >  < contig.fasta.fai >"
    exit(1) 

fs={}
fn={}
fl={}
contiglen={}

for l in open(sys.argv[2],"r"):
    ss=l.strip().split()
    contiglen[ss[0]]=int(ss[1])

for line in open(sys.argv[1],"r"):
    ss =line.strip().split()
    if len(ss) <=2 or line.startswith("#"):
        continue
    #numss = int(ss[2].replace("Contig",""))
    if ss[2] not in fs:
        fl[ss[2]] = contiglen[ss[2]]
        fs[ss[2]] = []
        fn[ss[2]] = 0
    fs[ss[2]].append(ss[0])
    fn[ss[2]] += 1 

cont = []
contigs=[]
sumn=0

for k,l in sorted(fl.items(), key=lambda item:item[1]):
    fv = 0 
    # add to list 
    for i in fs[k]:
        if i not in cont:
            cont.append(i)
            fv = 1
    if fv == 1 :
        contigs.append(k)
        sumn += fl[k]
        print "%s\t%s\t%s\t%s" %(k,sumn,len(cont),",".join(contigs))

    #    print "%s\t%s\t%s" %(k,contiglen[k],",".join(fs[k]))  
    #    sumn += contiglen[k]
#print "# Total: %s" %(sumn)
