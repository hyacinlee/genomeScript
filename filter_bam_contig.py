#!/usr/bin/python 
import sys

def main():
    keep=[]
    for line in open(sys.argv[1],"r"):
        keep.append(line.strip().split("\t")[0])

    for sam in sys.stdin:
        if sam.startswith("@"):
            if sam.startswith("@SQ") and sam.strip().split("\t")[1].replace("SN:","") not in keep:
                continue
            else:
                print sam.strip()
        else:
            if sam.strip().split("\t")[2] in keep:
                print sam.strip()

main()

