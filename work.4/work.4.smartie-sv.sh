#!/usr/bin/env bash
set -e
hostname
query=$1
target=$2
cpu=24
path=`pwd`
project=`echo $$`
export LD_LIBRARY_PATH=/work/xup/software/lib/htslib-1.10.2/lib:$LD_LIBRARY_PATH

### index
echo "Build target index"
cd ${path}
mkdir -p targetIndex && cd targetIndex
ln -fs ${target} target.fasta
/work/xup/software/Assemble/smartie-sv/bin/sawriter target.fasta
touch targetIndex.done

### align
echo "Aligning query to target"
cd ${path}
mkdir -p align && cd align
ln -s ${query} query.fasta
/work/xup/software/Assemble/smartie-sv/bin/blasr -clipping hard -alignContigs -sam -minMapQV 30 -nproc ${cpu} -minPctIdentity 50 -unaligned unaligned.fasta query.fasta ../targetIndex/target.fasta -out aligned.sam
### call
echo "Calling SVs"
cd ${path}
mkdir -p call && cd call
cat ../align/aligned.sam |/work/xup/software/Assemble/smartie-sv/bin/printgaps ../targetIndex/target.fasta all
echo "all done"
date
