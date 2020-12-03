#!/usr/bin/env python 
import sys 
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():

    args = get_args()

    if args.ft == "table":
        formatTable(args)

    if args.ft == "fasta":
        formatFasta(args)
      


def formatFasta(args):

    fasta = readFasta(args.f)
    for line in open(args.b,"r"):
        if line.startswith("#"):
            if args.add :
                print "%s\tSeq" % (line.strip()) 
            continue
        s         = line.strip().split()
        lociList  = joinPartList(s,args.bc).split("#")
        (fid,seq) = getSequence("".join(fasta[lociList[0]]),lociList)   
        if args.add :
            print "%s\t%s" % (line.strip(),seq)
        else: 
            print ">%s\n%s" % (fid,seq)





def formatTable(args):

    if not len(args.bc) == len(args.fc):
        print "Error, number of bait columns not eq the number of fish columns"
        exit(1)

    baitlist=[]
    for line in open(args.b,"r"):
        if line.startswith("#"):
            continue
        s=line.strip().split()
        tag=joinPartList(s,args.bc)
        baitlist.append(tag)

    #print baitlist
    for line in open(args.f,"r"):
        if line.startswith("#"):
            print line.strip()
            continue
        s=line.strip().split()
        tag=joinPartList(s,args.fc)
        #print tag
        if args.exc :
            if tag not in baitlist:
                print line.strip()
        else:
            if tag in baitlist:
                print line.strip()


def get_args():

    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton: use a table split file to handle fasta && table files 
    Writer  : Minghui Meng <hyacinlee@163.com> 

    1) Extract or Except saome line by mutiple-col number 
         formatTable.py -b bait.table -f fish.table
         formatTable.py -b bait.table -f fish.table --bc 1 2  --fc 2 3 --except  
    2) Add fasta seq to the end of table , 1 2 3 4 5 means chr|start|end|stand|name; 4 5 can be none
         formatTable.py -b bait.table -f fasta.file --bc 1 2 3 4 -add -ft fasta 
    3) Extract sub fasta ,
         formatTable.py -b bait.table -f fasta.file --bc 1 2 3 4 5 -ft fasta

'''
    )
    parser.add_argument('-b',metavar='',help='table split bait file ')
    parser.add_argument('-f',metavar='',help='fish file')
    parser.add_argument('-ft',metavar='',help='type of fish file  table|fasta ',choices=('table', 'fasta'),default="table")
    parser.add_argument('--bc',metavar='',help='set bait file columns,suppose like 1,2,3 ',nargs="+",default=[1])
    parser.add_argument('--fc',metavar='',help='set fish file columns,suppose like 2,3,4 ',nargs="+",default=[1])
    parser.add_argument('--exc',help='except bait lines',action='store_true',default=False)
    parser.add_argument('--add',help='add fasta to the end of table ',action='store_true',default=False)
    args = parser.parse_args()
    
    if not args.b or not args.f:
        parser.print_help()
        exit(1)

    return args




####    comman function 


def getSequence(seq,loci):
    mySeq = seq[int(loci[1])-1:int(loci[2])]
    myID  = "%s:%s-%s" % (loci[0],loci[1],loci[2])
    if len(loci) >= 4:
        if   loci[3] == "+":
            pass
        elif loci[3] == "-":
            mySeq = complementReverse(mySeq)
            myID  = "%s:%s-%s:%s" % (loci[0],loci[1],loci[2],loci[3])
        else:
            print "Error, %s is not +/- ,please choose the right columns number"
            exit(1)
    if len(loci) >= 5:
        myID = loci[4]

    return (myID,mySeq)



def joinPartList(mylist,colList):

    ss=[]
    for i in colList:
        ss.append(mylist[int(i)-1])

    return "#".join(ss)


def readTableDict(file,ck=0,cv=1):
    info={}
    with open(file,"r") as inf:
        for line in inf:
            ss = line.strip().split("\t")
            info[ss[ck]] = ss[cv]
    return info 


def readFasta(fasta_file):
    '''
        reading fasta file and return a dict of name and sequence .
    '''
    #print "# Reading fasta dict from %s" %(fasta_file)
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



def complementReverse(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.upper()
    return sequence[::-1]



def accumulateDict(myDict,val,keyA,keyB="no"):
    '''
        accumulate a 1~2d Dictionary  
    '''
    if keyB == "no":
        if keyA in myDict:
            myDict[keyA] += val
        else:
            myDict[keyA] = val
    else:
        if keyA in myDict:
            if keyB in myDict[keyA]:
                myDict[keyA][keyB] += val
            else:
                myDict[keyA].update({keyB: val})
        else:
            myDict.update({keyA:{keyB: val}})

    return myDict



if __name__ == '__main__':
    main()
