#!/usr/bin/env python 
from Bio import Entrez
import os,sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import sys, os, argparse, os.path,re,math,time
'''
database:
['pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest',
'structure', 'genome', 'books', 'cancerchromosomes', 'cdd', 'gap',
'domains', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pmc',
'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound',
'pcsubstance', 'snp', 'taxonomy', 'toolkit', 'unigene', 'unists']
'''

def main():

    args = getopt()
    dout=''
    if os.path.exists(args.out_dir):
        dout=os.path.abspath(args.out_dir)
    else:
        os.mkdir(args.out_dir)
        dout=os.path.abspath(args.out_dir)
    
    Entrez.email = "hyacinlee@163.com"     # Always tell NCBI who you are

    if args.term : 
        output_handle = open(dout+'/'+args.name+'.%s'%args.rettype, "w")
        handle = Entrez.esearch(db=args.database, term=args.term, idtype="acc")
        record = Entrez.read(handle)
        print record['IdList']
        print record
    
    if args.list:
        for line in open(args.list,"r"):
            i = line.strip().split()[0]

            handle = Entrez.efetch(db=args.database, id=i, rettype=args.rettype, retmode="text")
            record = SeqIO.read(handle, args.rettype)
            with open("%s/%s.fasta" %(dout,i),"w") as  out:
                SeqIO.write(record, out, args.rettype)

            print "Finish Fetch %s to file %s/%s.fasta" % (i,dout,i)




def getopt():
    parser = argparse.ArgumentParser(description='This script is used to fasta from ncbi ')
    parser.add_argument('-t','--term',help='input search  term : https://www.ncbi.nlm.nih.gov/books/NBK3837/#_EntrezHelp_Entrez_Searching_Options_',required=False)
    parser.add_argument('-l','--list',help='input list of ids ',required=False)
    parser.add_argument('-d','--database',help='Please input database to search nucleotide or protein  default nucleotide',default = 'nucleotide',required=False)
    parser.add_argument('-r','--rettype',help='return type fasta or gb default gb',default = "gb",required=False)
    parser.add_argument('-o','--out_dir',help='Please input  out_put directory path',default = os.getcwd(),required=False)
    parser.add_argument('-n','--name',default ='seq',required=False,help='Please specify the output, seq')
    args = parser.parse_args()
    
    if  not args.term and not args.list :
        exit(1)

    return args

if __name__ == '__main__':
    main()