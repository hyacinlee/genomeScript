#!/usr/bin/env python 
import sys 
import numpy as np
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass



def main():
    data={}
    order=[]
    args=get_args()
    rhflag=0
    for l in open(args.f,"r"):
        ls = l.strip().split("\t")
        if rhflag == 0 and args.rh:
            rhflag = 1 
            continue
        group=joinPartList(ls,args.b)   
        if group not in data:
             data[group]=[]
             order.append(group)
        data[group].append(float(ls[args.c-1]))
    print "#Group\tTotal\tCount\tMean\tMid\tMax\tMin"  
    for group in order:
         dl=data[group]
         print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %(group,sum(dl),len(dl),np.mean(dl),np.median(dl),max(dl),min(dl))
 



def joinPartList(mylist,colList):
    ss=[]
    for i in colList:
        if int(i) > len(mylist):
            ss.append("_")
            break
        else:
            ss.append(mylist[int(i)-1])

    return "#".join(ss)





def get_args():

    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton: cal mean, total , count , mid of dataset 
    Writer  : Minghui Meng <hyacinlee@163.com> 

         statTeblepy -b 1 2  -f table.file  -c 3  
'''
    )
    parser.add_argument('-b',metavar='',help='columns of group, link 1 2  ',nargs="+",default=[1])
    parser.add_argument('-f',metavar='',help='table file')
    parser.add_argument('-c',metavar='',help='columns of calcu',type=int,default=2)
    parser.add_argument('-rh',help='remove the first raw',action='store_true',default=False)
    args = parser.parse_args()
    
    if not args.b or not args.f:
        parser.print_help()
        exit(1)
    return args    



main()
