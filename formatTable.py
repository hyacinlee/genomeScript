#!/usr/bin/env python 
import sys 
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():

    args = get_args()

    if args.replace:
        replaceAll(args)
        exit(0)
    
    if args.order:
        order(args)
        exit(0)

    if args.f == "c" or args.f == "k" or args.f == "a" :
        otherFunction(args)
        exit(0)

    if args.ft == "table":
        formatTable(args)

    if args.ft == "fasta" and args.exc:
        removeFasta(args)

    if args.ft == "fasta" and args.link == -1 and not args.exc:
        formatFasta(args)

    if not args.link == -1 and args.ft == "fasta":
        linkfasta(args)


def otherFunction(args):
    count={}
    infos=[]
    adds={}

    for line in open(args.b):
        s  = line.strip().split("\t")
        ss = "\t".join([s[int(x)-1] for x in args.bc])
        if ss not in count:
            count[ss] = 1
            infos.append(line.strip())
            adds[ss]=[]
        else:
            count[ss] += 1

        if args.f == "a" :
            adds[ss].append(s[int(args.fc[0])-1])

    if args.f == "a" :
        for keys in count:
            print "%s\t%s" %(keys,",".join(adds[keys]))
    if args.f == "c" :
        for keys in count:
            print "%s\t%s" %(keys,count[keys])
    if args.f == "k" :
        print "\n".join(infos)


def order(args):
    fishDict = readTableDict(args.f,int(args.fc[0])-1,-1)
    for line in open(args.b):
        s  = line.strip().split("\t")
        ss = "\t".join([s[int(x)-1] for x in args.bc])
        if ss in fishDict:
            print fishDict[ss].strip()


def removeFasta(args):
    fasta,fo = readFasta(args.f,1)
    #print fo
    ex = []
    for line in open(args.b,"r"):
        s=line.strip().split("\t")[int(args.bc[0])-1]
        ex.append(s)
    for seqID in fo:
        if seqID in ex:
            continue
        else:
            print ">%s\n%s" %(seqID,"".join(fasta[seqID]))

            
def replaceAll(args):
    repDict={}
    for line in open(args.b,"r"):
        if line.startswith("#"):
            continue
        s=line.strip().split("\t")                
        repDict[s[int(args.bc[0])-1]] = s[int(args.bc[1])-1]

    if args.ft == "fasta":
        keep = 0 
        for line in open(args.f,"r"):
            if line.startswith(">"):
                fid = line.strip().replace(">","").split()[0]
                if fid in repDict:
                    print ">%s" % (repDict[fid])
                    keep = 1 
                else:
                    keep = 0
                    if not args.exc:
                        print ">%s" % (fid)
            else:
                if args.exc and keep == 0:
                    continue
                else:
                    #print ">%s" % (fid)
                    print line.strip()
    else:
        for line in open(args.f,"r"):
            if line.startswith("#"):
                print line.strip()
            else:
                s=line.strip().split("\t")
                repCol = s[int(args.fc[0])-1] 
                if repCol in repDict :
                    s[int(args.fc[0])-1] = repDict[repCol]
                    print "\t".join(s) 
                else:
                    if not args.exc :
                        print line.strip()


def linkfasta(args):

    loci={}     
    for line in open(args.b,"r"):
        if line.startswith("#"):
            continue
        else:
            s = line.strip().split("\t")
            qloci = [ s[int(i)-1] for i in args.bc ]
            if qloci[0] not in loci:
                loci[qloci[0]] = []
            #qlist = [ s[int(args.bc[1])-1],
            #qloci = [ s[int(i)-1] for i in args.bc ]
            loci[qloci[0]].append(qloci)

    #print loci
    fasta  = readFasta(args.f)
    splitN = "N"*args.link
    for chrid in loci.keys():
        seq = [ getSequence("".join(fasta[info[1]]),info[1:])[1]  for info in loci[chrid] ]
        print ">%s\n%s" %(chrid,"".join(seq))

#getSequence(seq,loci):



def formatFasta(args):

    fasta = readFasta(args.f)
 
    for line in open(args.b,"r"):
        if line.startswith("#"):
            if args.add :
                print "%s\tSeq" % (line.strip()) 
            continue
        s         = line.strip().split("\t")
        lociList  = joinPartList(s,args.bc).split("#")
        #print lociList[0] 
        #print type(fasta)
        #sys.stderr(lociList)
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
    baitinfo={}
    for line in open(args.b,"r"):
        if line.startswith("#"):
            continue
        s=line.strip().split("\t")
        tag=joinPartList(s,args.bc)
        baitlist.append(tag)
        baitinfo[tag] = line.strip()

    #print baitlist
    for line in open(args.f,"r"):
        if line.startswith("#"):
            print line.strip()
            continue
        s=line.strip().split("\t")
        #print s
        tag=joinPartList(s,args.fc)
        #print tag
        if args.exc :
            if tag not in baitlist:
                print line.strip()
        else:
            if tag in baitlist:
                if not args.add :
                    print line.strip()
                else:
                    print "%s\t%s" % (line.strip(),baitinfo[tag])
            else: # 
                if args.add:
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
         formatTable.py -b bait.table -f fish.table --bc 1 2  --fc 2 3 --exc
       Export both bait line and fish line
         formatTable.py -b bait.table -f fish.table --bc 1 2  --fc 2 3 --add 
    2) Add fasta seq to the end of table , 1 2 3 4 5 means chr|start|end|stand|name; 4 5 can be none
         formatTable.py -b bait.table -f fasta.file --bc 1 2 3 4 --add -ft fasta 
    3) Extract sub fasta ,
         formatTable.py -b bait.table -f fasta.file --bc 1 2 3 4 5 -ft fasta
    4) link fasta , 1 2 3 4 5 means NewChr| OldChr | OldStart | OldEnd | stand in NewChr 
         formatTable.py -b bait.table -f fasta.file -ft fasta --link 0 --bc 1 2 3 4 5
    5) replace sub colums or fasta ids
         formatTable.py -b bait.table -f fasta.file -ft fasta --bc 1 2  --replace  --exc
         formatTable.py -b bait.table -f fish.table -ft table --bc 3 2  --fc 4 --replace
    6) other function 
        a) count any colum of bait.table 
             formatTable.py -b tables  -f c --bc 1
        b) keep the first line of any colum 
             formatTable.py -b tables  -f k --bc 3 
    7) order fish file by bait file 
         formatTable.py -b tables -f fish.table --bc 3 --fc 2 --order 
'''
    )
    parser.add_argument('-b',metavar='xx',help='table split bait file ')
    parser.add_argument('-f',metavar='',help='fish file')
    parser.add_argument('-ft',metavar='',help='type of fish file  table|fasta ',choices=('table', 'fasta'),default="table")
    parser.add_argument('--bc',metavar='',help='set bait file columns,suppose like 1,2,3 ',nargs="+",default=[1])
    parser.add_argument('--fc',metavar='',help='set fish file columns,suppose like 2,3,4 ',nargs="+",default=[1])
    parser.add_argument('--exc',help='except bait lines',action='store_true',default=False)
    parser.add_argument('--order',help='order by  bait lines',action='store_true',default=False)
    parser.add_argument('--add',help='add fasta to the end of table ,or publish all bait info',action='store_true',default=False)
    parser.add_argument('--replace',help='replace table columns or fasta names',action='store_true',default=False)
    parser.add_argument('--link',help='link fasta with number of N',type=int,default=-1)
    args = parser.parse_args()
    
    if not args.b or not args.f:
        parser.print_help()
        exit(1)

    return args




####    common function 



def getSequence(seq,loci):
    #print loci
    if len(loci) == 1:
        return (loci[0],seq)

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
        if int(i) > len(mylist):
            ss.append("##")
            break
        else:
            ss.append(mylist[int(i)-1])

    return "#".join(ss)


def readTableDict(file,ck=0,cv=1):
    info={}
    with open(file,"r") as inf:
        for line in inf:
            ss = line.strip().split("\t")
            if cv != -1:
                info[ss[ck]] = ss[cv]
            else:
                info[ss[ck]] = line
    return info 


def readFasta(fasta_file,order=0):
    '''
        reading fasta file and return a dict of name and sequence .
    '''
    #print "# Reading fasta dict from %s" %(fasta_file)
    fo       = []
    fasta    = {}
    fasta_id = ''
    for line in open(fasta_file):
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
