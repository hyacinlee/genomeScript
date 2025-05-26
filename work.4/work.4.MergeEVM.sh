#!/usr/bin/sh 

# dir contian input.gff add.gff input.busco add.busco ref.fasta



grep mRNA input.gff |cut -f1,4,5,9 |sed 's/ID=//'|sed 's/;//'|sort -k1,1V -k2,2n > input.mRNA.bed
grep mRNA add.gff |cut -f1,4,5,9 |sed 's/ID=//'|sed 's/;//'|sort -k1,1V -k2,2n >   add.mRNA.bed
bedtools intersect -a input.mRNA.bed -b add.mRNA.bed -wo > Overlaps
awk '$2==$6 && $3==$7' Overlaps > Confirm.core.genes
awk '$2!=$6 || $3!=$7' Overlaps > Overlaps.candiate

formatTable.py -b Overlaps.candiate --bc 4 -f add.busco --fc 3 |wc -l
formatTable.py -b Overlaps.candiate --bc 8 -f add.busco --fc 3 |wc -l

formatTable.py -f Overlaps.candiate --fc 4 -b add.busco --bc 3 --add |perl -lnae 'if(@F==9){$F[3]="$F[3]\tnone\tnone"}else{$F[3]="$F[3]\t$F[9]\t$F[10]"};print join("\t",@F)' |cut -f1-10 > Overlaps.candiate.tmp1
formatTable.py -f Overlaps.candiate.tmp1 --fc 10 -b add.busco --bc 3 --add |perl -lane 'if(@F==10){print "$_\tnone\tnone"}else{print "$_"}' |cut -f1-12 > Overlaps.candiate.busco

awk '$5=="none" && $11!="none"' Overlaps.candiate.busco |cut -f1,2,3,4 |sort -u > Overlaps.candiate.busco.del
awk '$6=="Fragmented" && $11!="none"' Overlaps.candiate.busco |cut -f1,2,3,4 |sort -u >> Overlaps.candiate.busco.del

awk '$5=="none" && $11!="none"' Overlaps.candiate.busco |cut -f7-11    |sort -u > Overlaps.candiate.busco.tmp.add
awk '$6=="Complete" || $6=="Duplicated" ' Overlaps.candiate.busco > Overlaps.candiate.busco_for_add_to_del
formatTable.py -f Overlaps.candiate.busco.tmp.add -b Overlaps.candiate.busco_for_add_to_del --fc 4 --bc 10 --exc > Overlaps.candiate.busco.add


formatTable.py -b Overlaps --bc 8 -f add.mRNA.bed --fc 4 --exc > unone.Overlaps.bed
formatTable.py  -b unone.Overlaps.bed --bc 4 -f add.busco --fc 3 --add |perl -lnae 'print $_ if(@F>7)' |grep -v "#" |grep -v Fragmented|awk -v OFS="\t" '{print $6,$7,$8,$9,$1}' > unone.Overlaps.busco

cat unone.Overlaps.busco Overlaps.candiate.busco.add|sort -k1,1 -k2,2n > total.add.candiate
formatTable.py -f input.mRNA.bed --fc 4 -b Overlaps.candiate.busco.del --bc 4 --exc > input.mRNA.keep.bed
bedtools intersect -a input.mRNA.keep.bed -b total.add.candiate -wo > total.add.check 

perl -lne 'my @s=split(/=/,$_);$s[1]=~s/\;// ; print "$_\t$s[1]"' input.gff > input.gff.tmp
perl -lne 'my @s=split(/=/,$_);$s[1]=~s/\;// ; print "$_\t$s[1]"' add.gff > add.gff.tmp

formatTable.py -b input.mRNA.keep.bed --bc 4 -f input.gff.tmp --fc 10 |cut -f1-9 > result.keep.gff
formatTable.py -b total.add.candiate --bc 4 -f add.gff.tmp --fc 10    |cut -f1-9 |sed 's/;/_add;/'>    result.add.gff 

cat result.keep.gff result.add.gff |sort -V -k1,1 -k4,4n -k3,3r > result.fianl.gff
getGene.pl result.fianl.gff ref.fasta > result.final.cds.fasta 
cds2aa.pl result.final.cds.fasta  > result.final.pep.fasta

echo "/work/xup/bin/script/work.Busco.sh result.final.pep.fasta  /work/xup/database/BUSCO/actinopterygii_odb9/ prot"
