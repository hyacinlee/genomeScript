#!/usr/bin/sh 
gff=$1
awk '$3=="mRNA"' $1|sed 's/;/\t/g' | cut -f1,4,5,7,9 |sed ''s/ID=//|perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[4]\t0\t$F[3]"'
