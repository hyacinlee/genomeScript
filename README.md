# genomeScript

## gapFinder.py 
find gap in genome fasta 

## telBlock.py 
creat Telomere in chromsome start and end

## splitScaf.py 
split scaf and chromsome Fasta to contig by remove N

## mergeSV.sh
bash commond to merge sv result of  smartie and sniff

## formatTable.py 
usage: formatTable.py [-h] [-b] [-f] [-ft] [--bc  [...]] [--fc  [...]] [--exc]

RunCmd :
    Funciton: use a table split file to handle fasta && table files 
    Writer  : Minghui Meng <hyacinlee@163.com> 

    1) Extract or Except saome line by mutiple-col number 
        ./formatTable.py -b bait.table -f fish.table
        ./formatTable.py -b bait.table -f fish.table --bc 1,2  --fc 2,3 --except  

optional arguments:
  -h, --help    show this help message and exit
  -b            table split bait file (default: None)
  -f            fish file (default: None)
  -ft           type of fish file (default: table)
  --bc  [ ...]  set bait file columns,suppose like 1,2,3 (default: [1])
  --fc  [ ...]  set fish file columns,suppose like 2,3,4 (default: [1])
  --exc         except bait lines (default: False)
