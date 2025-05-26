#!/usr/bin/env python 
import sys
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():
    
    args=get_args()

    snplist=""
    if args.INPUT.endswith("vcf"):
        snplist=transVCF(args.INPUT)
    elif args.INPUT.endswith("lst"):
        snplist=args.INPUT
    else:
        print "\t#Error: Can't judge input type, input file must ends with .lst or .vcf "
        exit(1)
    
    (samplesInGroup,groupOfSample)=readGroup(args.g)

    out=open("result.groupAllele.%s_%s.xls" % (args.Q1,args.Q2),"w")

    sample=[]
    if  args.m == "none":
        out.write( "#Chr\tPos\tRef\tAlt\t%s\t%s\n" % (args.Q1,args.Q2) )
    else:
        out.write( "#Chr\tPos\tRef\tAlt\tGroup(%s)\tGroup(%s)\tGroup(%s)\t%s\n" % (args.Q1,args.Q2,args.m,"\t".join(samplesInGroup[args.m])))

    (t_number,q1_number,q2_number,ud_number)=(0,0,0,0)
    homAllen=["AA","GG","CC","TT"]

    for line in open(snplist,"r"):
        info=line.strip().split()

        if line.startswith("#"):
            sample = info[4:]
            #print sample
            checkGroup(sample,groupOfSample)
        else: 
            cid=info.pop(0)
            pos=info.pop(0)
            ref=info.pop(0)
            alt=info.pop(0)
            homAllen=[alt+alt,ref+ref]
            StatQ1 = statQ(info,sample,samplesInGroup[args.Q1])
            StatQ2 = statQ(info,sample,samplesInGroup[args.Q2])
            if args.hom:
                if StatQ1[0] not in homAllen or StatQ2[0] not in homAllen:
                    continue

            if "N" in StatQ1[0]  or "N" in StatQ2[0]:
                continue

            if StatQ1[0] == StatQ2[0] : 
                continue
            #print StatQ1,StatQ2,homAllen
            if float(StatQ1[1])/len(samplesInGroup[args.Q1]) < args.qr or float(StatQ2[1])/len(samplesInGroup[args.Q2]) < args.qr:
                continue
            
            t_number += 1
            if not args.m == "none":
                StatM = statQ(info,sample,samplesInGroup[args.m])
                #print StatQ1,StatQ2,homAllen,StatM
                if StatM[0] == "NN":
                    continue
                types="none"
                if StatM[0] ==StatQ1[0] and float(StatM[1])/len(samplesInGroup[args.m]) >= args.mr:
                    types=args.Q1
                    q1_number += 1
                elif StatM[0] ==StatQ2[0] and float(StatM[1])/len(samplesInGroup[args.m]) >= args.mr:
                    q2_number += 1
                    types=args.Q2
                else:
                    ud_number += 1
                    types = "Undef"

                mdQ = (getGene(info,sample,samplesInGroup[args.m]))[1]
                mGene = map(lambda x:mdQ[x],samplesInGroup[args.m])    

                out.write( "%s\t%s\t%s\t%s\t%s(%s/%s)\t%s(%s/%s)\t%s(%s/%s)\t%s\n" % (cid,pos,ref,alt,StatQ1[0],StatQ1[1],len(samplesInGroup[args.Q1]),StatQ2[0],StatQ2[1],len(samplesInGroup[args.Q2]),types,StatM[1],len(samplesInGroup[args.m]),"\t".join(mGene)))
            else:
                out.write(  "%s\t%s\t%s\t%s\t%s(%s/%s)\t%s(%s/%s)\n" % (cid,pos,ref,alt,StatQ1[0],StatQ1[1],len(samplesInGroup[args.Q1]),StatQ2[0],StatQ2[1],len(samplesInGroup[args.Q2])))


    out.close()

    if not args.m == "none":
        print  "# Result:\n# Total population-specific SNPs: %s \n# Mid-pop same with Q1: %s\n# Mid-pop same with Q2: %s\n# Mid-pop undefine: %s" % (t_number,q1_number,q2_number,ud_number)



def transVCF(file):

    snplist="./result.tmp.trans.lst"
    out=open(snplist,"w")

    for line in open(file,"r"):
        if line.startswith("##"):
            continue

        info=line.strip().split() 
        if line.startswith("#CHROM"):       
            out.write("%s\t%s\t%s\t%s\t%s\n" % (info[0],info[1],info[3],info[4],"\t".join(info[9:])))
        else:
            gts=[]
            for ss in info[9:]:
                gt=ss.split(":")[0]
                if gt == "0/0" or gt == "0|0":
                    gts.append("%s%s"%(info[3],info[3]))
                elif gt == "0/1"or gt == "0|1" or gt == "1|0":    
                    gts.append("%s%s"%(info[3],info[4]))
                elif gt == "1/1" or gt == "1|1":    
                    gts.append("%s%s"%(info[4],info[4]))
                elif gt == "./." or gt == ".|.":    
                    gts.append("NN")
                else:
                    print "# %s is not a right format of GenoType in vcf file, at %s %s" %(gt,info[0],info[1])

            out.write("%s\t%s\t%s\t%s\t%s\n" % (info[0],info[1],info[3],info[4],"\t".join(gts)))

    return snplist 


def statQ(info,sample,Q):
    # return the max allen
    stat=["none","none"]
    count={}
    Q=(getGene(info,sample,Q))[0]
    stat=judgeMin(Q)  
    return stat


def getGene(info,sample,Q):

    gQ=[]
    dgQ={}
    for i in range(len(info)):
        if sample[i] in Q:
            gQ.append(info[i])
            dgQ[sample[i]] = info[i]
    return gQ,dgQ


def  checkGroup(sample,groupOfSample):

    for sampleName in groupOfSample.keys():
        if sampleName not in sample:
            print "#Error! Sample ID %s in group file not in snplist file !" % (sampleName)
            exit(1) 


def readGroup(file):
    group={}
    sample={}

    for line in open(file,"r"):
        if line.startswith("#"):
            continue
        info=line.strip().split()
        if info[1] not in group:
            group[info[1]]=[]

        group[info[1]].append(info[0])
        sample[info[0]] = info[1]

    return (group,sample)



def accumulateDict(myDict,val,keyA,keyB="no"):
    '''
        accumulate a 1~2d Dictionary  
    '''
    if keyB == "no":
        if keyA in myDict:
            myDict[keyA] += val
        else:
            myDict[keyA] = val
    else:
        if keyA in myDict:
            if keyB in myDict[keyA]:
                myDict[keyA][keyB] += val
            else:
                myDict[keyA].update({keyB: val})
        else:
            myDict.update({keyA:{keyB: val}})

    return myDict


def judgeMin(mylist):
    count={}
    for ss in mylist:
        if ss not in count:
            count[ss]=1
        else:
            count[ss]+=1
    a = sorted(count.items(), key=lambda x: x[1],reverse=True)
    return a[0]


def get_args():

    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton:  Identification of  population-specific SNPs.
    Writer:    Meng Minghui < hyacinlee@163.com >
'''
    )
    parser.add_argument('-g',metavar='group',help='input of group file',type=str)
    parser.add_argument('-Q1',metavar='Q1',help='name of group1',type=str)
    parser.add_argument('-Q2',metavar='Q2',help='name of group2',type=str)
    parser.add_argument('-m',metavar='m',help='ame of middle group , can be false',default="none")
    parser.add_argument('-qr',metavar='qr',help='min-rate to define the SNP type of Q1 and Q2 ',default="0.8",type=float)
    parser.add_argument('-mr',metavar='mr',help='min-rate to define the SNP type of Mid group, use with -m ',default="0.7",type=float)
    parser.add_argument('-hom',help='only keep hom allen for parents',action='store_true',default=False)
    parser.add_argument('INPUT',metavar='input',help='cmd input vcf file must be given')


    args = parser.parse_args()
    return args



if __name__ == '__main__':
    main()
