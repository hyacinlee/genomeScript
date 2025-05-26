#!/usr/bin/env python

import sys

if ('-h' in sys.argv) or ('--help' in sys.argv) or len(sys.argv) < 2:
    print >> sys.stderr, sys.argv[0] + ' works for change ID/Parent in gff\n' + \
        sys.argv[0] + ' id.o2n.list[old:1, new:2] old.gff > new.gff'
    exit(0)

ids = {}

with open(sys.argv[1], 'r') as fI:
    for line in fI:
        idO, idN = line.strip().split()
        ids[idO] = idN

with open(sys.argv[2], 'r') as fG:
    for line in fG:
        infos = line.strip().split()
        out = infos[:8]
        out_d = ''
        for i in infos[8].strip(';').split(';'):
            items = i.split('=')
            if items[0] == 'ID':
                out_d += 'ID=' + ids[items[1]] + ';'
            elif items[0] == 'Parent':
                if items[1] not in ids: continue
                out_d += 'Parent=' + ids[items[1]] + ';'
            elif items[0] != 'Name':
                out_d += i + ';'
        out.append(out_d)
        print >> sys.stdout, '\t'.join(out)
