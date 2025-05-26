#!/usr/bin/sh 

# dir contian input.gff add.gff input.busco add.busco ref.fasta

grep mRNA input.gff |cut -f1,4,5,9 |sed 's/ID=//'|sed 's/;//'|sort -k1,1V -k2,2n > input.mRNA.bed
grep mRNA add.gff |cut -f1,4,5,9 |sed 's/ID=//'|sed 's/;//'|sort -k1,1V -k2,2n >   add.mRNA.bed
#bedtools intersect -a input.mRNA.bed -b add.mRNA.bed -wo > Overlaps

### keep core busco of input 
awk '$2=="Complete" || $2=="Duplicated" ' input.busco > input.busco.kept 
echo "# Orgin busco number:  `cut -f1 input.busco.kept |sort -u |wc -l`"
awk '$2=="Complete" || $2=="Duplicated" ' add.busco > add.busco.raw

formatTable.py -b input.busco.kept -f add.busco.raw --exc > add.busco.candidated
formatTable.py -b add.busco.candidated -f add.mRNA.bed --bc 3 --fc 4 >  add.busco.candidated.genes 
formatTable.py -b add.busco.candidated -f add.mRNA.bed --bc 3 --fc 4 --add |perl -lnae 'print "$_" if(@F>5)' > add.busco.candidated.genes.info 
echo "# Add busco in gff2: `cut -f1 add.busco.candidated|sort -u |wc -l`"

# check overlap 
bedtools intersect -a input.mRNA.bed -b add.busco.candidated.genes -wo > Overlaps 
formatTable.py -b Overlaps -f add.busco.candidated.genes.info --bc 8 --fc 4 --exc > Overlaps.noneOver.adds 
cut -f4 Overlaps.noneOver.adds > Overlaps.noneOver.adds.ids 
echo "# Add busco in gff2 unOverlaps: `cut -f5  Overlaps.noneOver.adds|sort -u |wc -l`"


formatTable.py -f Overlaps --fc 4 -b input.busco --bc 3 --add |perl -lnae 'if(@F==9){$F[3]="$F[3]\tnone\tnone"}else{$F[3]="$F[3]\t$F[9]\t$F[10]"};print join("\t",@F)'|cut -f1-10>Overlaps.tmp1  
formatTable.py -f Overlaps.tmp1 --fc 10 -b add.busco --bc 3 --add |perl -lane 'if(@F==10){print "$_\tnone\tnone"}else{print "$_"}' |cut -f1-12 > Overlaps.raw

awk '$6=="none" || $6=="Fragmented"' Overlaps.raw > Overlaps.candiate
cut -f4 Overlaps.candiate |sort -u > Overlaps.input.dels 
cut -f10 Overlaps.candiate |sort -u > Overlaps.adds
formatTable.py -b Overlaps.adds -f add.busco --bc 1 --fc 3|grep -v "#" > Overlaps.adds.busco
echo "# Add busco in gff2 Overlaps: `cut -f1  Overlaps.adds.busco|sort -u |wc -l`"

cat Overlaps.noneOver.adds.ids Overlaps.adds > Overlaps.final.add.ids 

 
#formatTable.py -b Overlaps -f add.busco.candidated.genes --bc 8 --fc 4 --exc >> Overlaps.adds 

fishInWinter.pl -bf table -ff gff -except Overlaps.input.dels input.gff > input.del.gff
fishInWinter.pl -bf table -ff gff Overlaps.final.add.ids add.gff > add.tmp.gff 
grep mRNA add.tmp.gff |cut -f9 |uniq |sed 's/ID=//' |sed 's/;//'|perl -lnae '$s++;$num=(sprintf "%05d",$s);print "$F[0]\tADD$num"' > add.tmp.rename
/work/xup/bin/script/gffIdChange.py add.tmp.rename add.tmp.gff > add.add.gff

cat add.add.gff input.del.gff |sort -k1,1V -k4,4n -k3,3r > result.tmp.gff 
grep mRNA  result.tmp.gff |cut -f9 |uniq |sed 's/ID=//' |sed 's/;//'|perl -lnae '$s++;$num=(sprintf "%05d",$s);print "$F[0]\tMERGE00$num"' > result.mRNA.change
/work/xup/bin/script/gffIdChange.py result.mRNA.change result.tmp.gff > result.final.gff

getGene.pl result.final.gff ref.fasta > result.final.cds.fasta
cds2aa.pl result.final.cds.fasta  > result.final.pep.fasta

echo "/work/xup/bin/script/work.Busco.sh result.final.pep.fasta  /work/xup/database/BUSCO/actinopterygii_odb10/ prot" > busco.sh 

