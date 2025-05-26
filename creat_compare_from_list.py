#!/usr/bin/env python 
import sys 

lst=[]

for l in open(sys.argv[1],"r"):
    lst.append(l.strip().split()[0])
    

compared=[]
out1=open("Compare.lst","w")
out2=open("Compare.all.lst","w")
for x in lst:
    out2.write(x+"\t"+x+"\n")
    for y in lst:
        if [x,y] not in compared and not x == y :
            out1.write(x+"\t"+y+"\n")
            out2.write(x+"\t"+y+"\n")
            out2.write(y+"\t"+x+"\n")
            compared.append([x,y])
            compared.append([y,x])

out1.close()
out2.close()
