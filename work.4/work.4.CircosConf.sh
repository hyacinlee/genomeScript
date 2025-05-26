#!/usr/bin/sh 

if [ $# -lt 1 ]
then
        echo "Usage:circos.conf.sh  < input dir >  < windows size>";
        echo "input file:"
        echo -e "\tGenome.fasta  Gene.gff  Repeat.gff  | *bams | sv.bed | snp.bed | indel.bed | coords(link)"
        exit;
fi

#numcer="/work/xup/software/Base/mummer-4.0.0beta2/bin/"
BIN=$(cd `dirname $0`;pwd )


input=$1
windows=$2

if [ ! -e circos.conf ];then 
   cp -r $BIN/circos/* .
fi 

mkdir -p data && cd data  

ln -sf $input/Genome.fasta &&  samtools faidx Genome.fasta

perl $BIN/get_GC_by_windows.pl  $input/Genome.fasta $windows GC.txt
grep Chr Genome.fasta.fai |sort -k1V |cut -f1,2 |perl -lane 'print "chr\t-\t$F[0]\t$F[0]\t0\t$F[1]\tlgrey"' > karyotype.txt 
grep Chr Genome.fasta.fai |sort -k1V |perl -lane 'print "$F[0]\t1\t$F[1]"'> chr.bed
#cat Genome.fasta.fai |cut -f1,2 |perl -lane 'print "chr\t-\t$F[0]\t$F[0]\t0\t$F[1]\tlgrey"' > karyotype.txt
#cat Genome.fasta.fai |perl -lane 'print "$F[0]\t1\t$F[1]"'> chr.bed
bedtools makewindows -b chr.bed -w $windows > windows.bed

grep mRNA $input/Gene.gff |cut -f1,4,5 > gene.bed 
bedtools coverage -a windows.bed -b gene.bed |cut -f1,2,3,7 > gene.density

grep Gypsy $input/Repeat.gff | cut -f1,4,5 > Gypsy.bed 
grep Copia $input/Repeat.gff | cut -f1,4,5 > Copia.bed 
cat $input/Repeat.gff | cut -f1,4,5 > repeat.bed
bedtools coverage -a windows.bed -b repeat.bed |cut -f1,2,3,7 > repeat.density
bedtools coverage -a windows.bed -b Gypsy.bed  |cut -f1,2,3,7 > Gypsy.density 
bedtools coverage -a windows.bed -b Copia.bed  |cut -f1,2,3,7 > Copia.density

exit

ls  $input/*bam |while read i;do 
    na=`basename $i .bam`
    echo "bedtools coverage -mean -b ${i} -a windows.bed > ${na}.depth"
done > cmd.depth.sh 

perl -lane 'print "$F[0]\t$F[1]\t$F[1]"' $input/snp.info |grep -v "#" > snp.bed
bedtools coverage -a windows.bed -b snp.bed | cut -f1,2,3,4 > snp.density
perl -lane 'print "$F[0]\t$F[1]\t$F[1]"' $input/indel.info |grep -v "#" > indel.bed
bedtools coverage -a windows.bed -b indel.bed | cut -f1,2,3,4 > indel.density
bedtools coverage -a windows.bed -b sv.bed | cut -f1,2,3,4 > sv.density

awk '$7>30000 && $10>90 && $15!=$16' $input/coords |perl -lane 'my $s=$F[14];$s=~s/Chr//;$s=~s/^0//;print "$F[14]\t$F[0]\t$F[1]\t$F[15]\t$F[3]\t$F[4]\tcolor=chr$s"' > link.txt

