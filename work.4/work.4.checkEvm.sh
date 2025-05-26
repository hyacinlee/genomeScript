#!/usr/bin/sh 
set -e 
gff=$1
cds=$2
name=`basename $gff .gff`

#getGene.pl $gff $genome > ${name}.cds
cds2aa.pl --check $cds > ${cds}.check
cat ${cds}.check  |sed '1d' |perl -lane 'if($F[3]*$F[4]==0){print $F[0]}' > ${name}.uncorrect.ids
fishInWinter.pl -bf table -ff gff --except ${name}.uncorrect.ids $gff   > ${name}.check.gff 
fishInWinter.pl -bf table -ff fasta --except ${name}.uncorrect.ids $cds > ${name}.check.cds.fasta 
cds2aa.pl ${name}.check.cds.fasta > ${name}.check.pep.fasta
#fishInWinter.pl -bf table -ff fasta --except ${1}.uncorrect.ids $1.evm.tran.pep > ${name}.check.pep.fasta
echo "/work/xup/bin/src/work.Busco.sh ${name}.check.pep.fasta /work/xup/database/BUSCO/actinopterygii_odb10/ prot" # busco.sh 
#qsub -cwd -pe smp 24 busco.sh 
