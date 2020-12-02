#cat ../sniffles.vcf | perl -F"\t" -lane'if(/END=(.*?);/){$end=$1+1;if(/SVTYPE=(.*?);/){$sv=$1; next if($sv eq  "BND") ;if (/SVLEN=(.*?);/){$len=abs($1);print "$F[0]\t$F[1]\t$end\t$len\t$sv"}}}' |sed 's/Chromosome/GRCz11_/' > sniffles.sv.bed


#bedtools intersect -a sniffles.sv.bed -b all.svs.bed -wo  > merge.bed


#perl -lane 'if(/#/){next}; my $l=length($F[4]);if($F[3] eq "N" && $l >= 50){my $s=my $e=$F[1]+$l-1; print "$F[0]\t$F[1]\t$e\t$l\tINS"}' ../sniffles.vcf


perl -lane 'if(/#/){next}; my $l=length($F[4]);if($F[3] eq "N" && $l >= 50){my $s=$F[1]-$l;my $e=$F[1]+$l; print "$F[0]\t$s\t$e\t$F[1]\t$l\tINS"}' ../sniffles.vcf|sed 's/Chromosome/GRCz11_/'   > sniffles.INS.bed
perl -lane 'if(/END=(.*?);/){$end=$1;if(/SVTYPE=(.*?);/){$sv=$1; next if($sv eq  "BND") ;if (/SVLEN=(.*?);/){$len=abs($1);print "$F[0]\t$F[1]\t$end\t$len\t$sv"}}}' ../sniffles.vcf |grep DEL | sed 's/Chromosome/GRCz11_/' > sniffles.DEL.bed

grep insertion ../all.svs.bed | sed 's/Chromosome/GRCz11_/'   > samartie.INS.bed
grep del ../all.svs.bed       | sed 's/Chromosome/GRCz11_/'   > samartie.DEL.bed



bedtools intersect -a sniffles.INS.bed -b samartie.INS.bed -wo > tmp.merge.INS.bed 
cut -f1,2,3 tmp.merge.INS.bed  |uniq -c |awk '$1!=1' | perl -lane 'print "$F[1]\t$F[2]\t$F[3]"' > tmp.INS.remove1.ids
cut -f7,8,9 tmp.merge.INS.bed | uniq -c |awk '$1!=1' | perl -lane 'print "$F[1]\t$F[2]\t$F[3]"' > tmp.INS.remove2.ids
./formatTable.py -b tmp.INS.remove1.ids -f tmp.merge.INS.bed       --bc 1 2 3 --fc 1 2 3 --exc    > tmp.INS.remove1.bed
./formatTable.py -b tmp.INS.remove2.ids -f tmp.INS.remove1.bed --bc 1 2 3 --fc 7 8 9 --exc    > tmp.INS.remove2.bed
echo -e "#chr1\tstart1\tend1\tType\tlength\tchr2\tstart2\tend2" > merge.filter.INS.bed 
perl -lane 'print "$F[0]\t$F[3]\t$F[3]\t$F[5]\t$F[10]\t$F[12]\t$F[13]\t$F[14]"' tmp.INS.remove2.bed >> merge.filter.INS.bed


bedtools intersect -a sniffles.DEL.bed -b samartie.DEL.bed -wo > tmp.merge.DEL.bed 
cut -f1,2,3 tmp.merge.DEL.bed  |uniq -c |awk '$1!=1' | perl -lane 'print "$F[1]\t$F[2]\t$F[3]"' > tmp.DEL.remove1.ids
cut -f6,7,8 tmp.merge.DEL.bed | uniq -c |awk '$1!=1' | perl -lane 'print "$F[1]\t$F[2]\t$F[3]"' > tmp.DEL.remove2.ids
./formatTable.py -b tmp.DEL.remove1.ids -f tmp.merge.DEL.bed       --bc 1 2 3 --fc 1 2 3 --exc    > tmp.DEL.remove1.bed
./formatTable.py -b tmp.DEL.remove2.ids -f tmp.DEL.remove1.bed --bc 1 2 3 --fc 6 7 8 --exc    > tmp.DEL.remove2.bed
echo -e "#chr1\tstart1\tend1\tType\tlength\tchr2\tstart2\tend2" > merge.filter.DEL.bed 
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[3]\t$F[11]\t$F[12]\t$F[13]"' tmp.DEL.remove2.bed >> merge.filter.DEL.bed
