#!/usr/bin/sh 
compare=$1
indir=$2

# indir:
# *.cds    : cds.fasta 
# *.seqids : chromosome ids 
# *.bed    : gene bed file ,use gff2bed.sh *gff > *bed 

cat $compare | sed 's/\t/\n/g' | uniq > input.lst

cat input.lst |while read i;do ln -sf ${indir}/${i}.bed ${i}.bed;done
cat input.lst |while read i;do ln -sf ${indir}/${i}.cds ${i}.cds ;done

#rename .Chr  Chr *bed *cds

#cat compare.txt | sed 's/\.//g' > compare.rename.txt

cat $compare |while read a b ;do echo "python -m jcvi.compara.catalog ortholog ${a} ${b}";done  | parallel -j 8

cat $compare |while read a b ;do echo "python -m jcvi.compara.synteny screen --minspan=25 --simple ${a}.${b}.anchors ${a}.${b}.anchors.new";done | parallel -j 8

cat input.lst |while read i;do cut -f1 ${indir}/${i}.seqids |paste -sd, ;done  > seqids

cat input.lst |perl -lnae 'BEGIN{$s=1} print ".$s,.15,.85,0,m,$F[0],top,$F[0].bed";$s++' > layout
cat $compare  |perl -lnae 'BEGIN{$s=0} my $e=$s+1;print "e,$s,$e,$F[0].$F[1].anchors.simple";$s++' >> layout


echo "python -m jcvi.graphics.karyotype seqids layout"

#cat input.lst |while read i;do 

#python -m jcvi.compara.synteny screen --minspan=30 --simple bra.ath.anchors bra.ath.anchors.new
#python -m jcvi.compara.synteny screen --minspan=30 --simple dbc.bra.anchors dbc.bra.anchors.new
#python -m jcvi.compara.synteny screen --minspan=30 --simple ath.dbc.anchors ath.dbc.anchors.new


#rm -rf *des *prj *sds *ssp *suf *tis
