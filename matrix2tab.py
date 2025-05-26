#!/usr/bin/env python
import sys 


head=[]
for line in open(sys.argv[1],"r"):
    s = line.strip().split()
    if len(head) == 0 :
        head = s
    else:
        sam = s.pop(0)
        for x in range(len(s)):
            print "%s\t%s\t%s" %(sam,head[x],s[x])    
 
