#################################################################
# take a dirtory of sam files and call all cigar code variants
# in each sam file

## author: Andrew Ruttenberg
## contact: ruttenberg.andrew@wustl.edu
#################################################################

import argparse
import os
from collections import Counter

class indel:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom=chrom
        self.pos=pos
        self.alt=alt
        self.ref=ref
    
    def __str__(self):
        return(f"{self.chrom}:{self.pos}  {self.ref}  {self.alt}")

    def __eq__(self, other):
        return self.chrom==other.chrom and self.pos==other.pos and self.alt==other.alt and self.ref==other.ref

class indelCounter:
    def __init__(self, indel, count):
        self.indel=indel
        self.count=count

    def __eq__(self, other):
        return self.indel==other.indel

    def __hash__(self):
        return hash((self.indel.chrom, self.indel.pos))
        

def makeFastaDict(fasta):
    with open(fasta, 'r') as fa:
        lines = [l for l in fa]
    firstLine=lines.pop(0)
    chrom=firstLine.split()[0][1:]
    sequence=""
    FastaDic = {}
    for line in lines:
        if line.startswith(">"):
            FastaDic[chrom]=sequence
            chrom=line.split()[0][1:]
            sequence=""
        else:
            sequence=sequence+line.replace('\n', '')
    FastaDic[chrom]=sequence
    return(FastaDic)

def cigarToList(cigarCode):
    out=[]
    workingNumber=""
    for i in range(len(cigarCode)):
        if cigarCode[i].isdigit():
            workingNumber=workingNumber+cigarCode[i]
        else:
            cigType=cigarCode[i]
            if cigType=='=' or cigType=='X':
                cigType='M'    
            out.append([cigType, int(workingNumber)])
            workingNumber=""
    return(out)

def makeVCF(indels, VCFPath):
    f=open(VCFPath, "w")
    with open("/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/Indel/CigarOutput/References/header.txt", "r") as header:
        for line in header:
            f.write(line)
    f.write(f'##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of variant"\n')
    f.write(f'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write(f'##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
    f.write(f'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    for key in indels:
        FinalIndel=indels[key]
        f.write(f"{FinalIndel.indel.chrom}\t{FinalIndel.indel.pos}\t.\t{FinalIndel.indel.ref}\t{FinalIndel.indel.alt}\t.\t.\tCount={FinalIndel.count};TYPE=INDEL\tGT:GQ:DP\t1/1:.:.\n")
    f.close()

def detectIndels(samFile, ref, currentIndelDict, LRS, indels):
    indelDict=dict()

    with open(samFile, "r") as sam_file:
        for line in sam_file:
            info=line.strip().split('\t')
            cigar=info[5]
            if "I" in cigar or "D" in cigar:
                chrom=info[2]
                pos=int(info[3])-1
                seq=info[9]
                spot=-1
                cigarList=cigarToList(cigar)
                for i in range(len(cigarList)):
                    el=cigarList[i]
                    if el[0]=='I' and i==0:
                        REF=ref[pos-1].replace('R', 'N')
                        ALT=ref[pos-1].replace('R', 'N')+seq[spot+1:spot+el[1]+1].replace('R', 'N')
                        if not indels and len(ALT)<=50:
                            spot=spot+el[1]
                            continue
                        if indels and len(ALT)>50:
                            spot=spot+el[1]
                            continue
                        newIndel=indelCounter(indel(chrom, pos, REF, ALT), 1)
                        if str(newIndel.indel) in indelDict:
                            indelDict[str(newIndel.indel)].count+=1
                        else:
                           indelDict[str(newIndel.indel)]=newIndel
                        spot=spot+el[1]
                        continue
                    elif el[0]=='M':
                        pos=pos+el[1]
                        spot=spot+el[1]
                    elif el[0]=="S" or el[0]=="H":
                        break
                    elif el[0]=="I":
                        REF=ref[pos-1].replace('R', 'N')
                        ALT=ref[pos-1].replace('R', 'N')+seq[spot+1:spot+el[1]+1].replace('R', 'N')
                        if not indels and len(ALT)<=50:
                            spot=spot+el[1]
                            continue
                        if indels and len(ALT)>50:
                            spot=spot+el[1]
                            continue
                        newIndel=indelCounter(indel(chrom, pos, REF, ALT), 1)
                        if str(newIndel.indel) in indelDict:
                            indelDict[str(newIndel.indel)].count+=1
                        else:
                           indelDict[str(newIndel.indel)]=newIndel
                        spot=spot+el[1]
                    elif el[0]=="D":
                        REF=ref[pos-1:pos+el[1]].replace('R', 'N')
                        ALT=seq[spot].replace('R', 'N')
                        if not indels and len(REF)<=50:
                            pos=pos+el[1]
                            continue
                        if indels and len(REF)>50:
                            pos=pos+el[1]
                            continue
                        newIndel=indelCounter(indel(chrom, pos, REF, ALT), 1)
                        if str(newIndel.indel) in indelDict:
                            indelDict[str(newIndel.indel)].count+=1
                        else:
                           indelDict[str(newIndel.indel)]=newIndel
                        pos=pos+el[1]
                    else:
                        break
        for key in indelDict:
            if key==currentIndelDict:
                currentIndelDict[key].count=max(currentIndelDict[key].count, indelDict[key].count)    
            else:
                currentIndelDict[key]=indelDict[key]
        return(currentIndelDict) 
def main(): 

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--samDir', required=True, help='inputed sam files for each window')
    parser.add_argument('-r', '--ref', required=True, help='Ref seq')
    parser.add_argument('-o', '--outfile', required=True, help='name of outfile')
    parser.add_argument('-l', '--LRS', required=True, help='bool for if the input is LRS or not')
    parser.add_argument('-i', '--Indel', required=True, help='bool for if detecting indels or SVs (true means indels)')
    args = parser.parse_args()
    print("making ref dictonary")
    FastaDict=makeFastaDict(args.ref)
    TotalIndelList=dict()
    SamFiles=os.listdir(args.samDir)
    numFile=len(SamFiles)
    i=0
    print("running pileup")
    for filename in os.listdir(args.samDir):
        i+=1
        if i % 1000 == 0:
            print(f"{i} files read so far")
        file_path = os.path.join(args.samDir, filename)
        if os.path.isfile(file_path):
            basename=filename.split('.')[0]
            chrom=basename.split('_')[0]
            refSeq=FastaDict[chrom]
            temp=detectIndels(file_path, refSeq, TotalIndelList, eval(args.LRS), eval(args.Indel))
            TotalIndelList=temp
    makeVCF(TotalIndelList, args.outfile)


if __name__ == '__main__':
    main()
