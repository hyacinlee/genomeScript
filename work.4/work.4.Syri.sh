#!/usr/bin/sh 
set -e 
ref="./JT.chromosome.fasta"
query="./Dongyajiangtun.fasta"
name="Dongyajiangtun"

#formatTable.py -b Dongya.idchange -f Dongyajiangtun.fasta.fai --bc 1 --fc 1 --add |perl -lnae 'print "$F[0]\t1\t$F[1]\t$F[7]\t$F[6]"' > Dongyajiangtun.idchange.bed
#formatTable.py -b Dongyajiangtun.idchange.bed -f Dongyajiangtun.fasta --bc 1 2 3 4 5 -ft fasta > Dongyajiangtun.idchange.fasta
/work/xup/software/Assemble/minimap2-master/minimap2 $ref $query --eqx -a -x asm5 -o ${name}.sam -t 24
#/work/xup/software/Base/mummer-4.0.0beta2/bin/show-coords -THrd Dongyajiangtun.delta > Dongyajiangtun.coords
unset PYTHONPATH
source activate py36
python3 /work/xup/software/Pan/syri-1.4/syri/bin/syri -c ${name}.sam -r $ref -q $query -k -F S 
python3 /work/xup/software/Pan/syri-1.4/syri/bin/plotsr syri.out $ref $query -H 10 -W 10
