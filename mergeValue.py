#!/usr/bin/env python 
import sys 

def main():
    if len(sys.argv) < 3:
        #print sys.argv
        print "\tUsage: %s < in tables > < sample cols, can be mutiple > "
        print "\tin tables:  [Name]  [value col]  [file name]"
        print sys.argv
        
        exit(1)


    info={}
    fs=[]

    for l in open(sys.argv[1],"r"):
        [name,col,files] =  l.strip().split()
        fs.append(name)
        for line in open(files,"r"):
            ls = line.strip().split()
            va = ls[int(col)-1]
            k = joinPartList(ls,sys.argv[2:])
            if k not in info:
                info[k]={}
            info[k][name] = va 

    print "Info\t%s" %("\t".join(fs))
    for k in info.keys():
        ol = []
        for name in fs:
            if name not in info[k]:
                info[k][name]="NA"
            ol.append(info[k][name])
        print "%s\t%s" % (k,"\t".join(ol))


def joinPartList(mylist,colList):
    ss=[]
    for i in colList:
        if int(i) > len(mylist):
            ss.append("##")
            break
        else:
            ss.append(mylist[int(i)-1])

    return "#".join(ss)

main()
