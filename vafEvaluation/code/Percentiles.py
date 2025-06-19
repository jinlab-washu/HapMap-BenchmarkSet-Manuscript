import argparse
import matplotlib.pyplot as plt
import numpy as np

class Allele:
    def __init__ (self, base, count, depth):
        self.base=base
        self.count=count
        self.depth=depth

    def __eq__ (self, other):
        return self.base==other.base
    def __str__ (self):
        return f"base: {self.base}, freq: {self.count/self.depth}"

class SNV:
    def __init__ (self, chrom, POS, REF, ALT, line):
        self.chrom=chrom
        self.POS=POS
        self.REF=REF
        self.ALT=ALT
        self.line=line


def makeSNVs(truthset, Calls):
    OutDict={}

    with open(truthset, 'r') as truthfile:
        linesTruth = [l for l in truthfile if not l.startswith('#') and not l.startswith('"##')]
    for line in linesTruth:
        split=line.split()
        chrom=split[0].split('r')[1]
        if chrom=='M':
            continue
        if chrom=='X':
            chrom=23
        elif chrom=='Y':
            chrom=24
        else:
            chrom=int(chrom)
        POS=int(split[1])
        REF=split[3]
        ALTBase=split[4]
        AFString=split[9].split(':')
        if len(AFString)==1:
            AF=0
        elif AFString[1]=='.':
            AF=0
        else:
            AF=float(AFString[1])
        ALT=Allele(ALTBase, AF, 1)
        key=f"{chrom}_{POS}_{REF}_{ALT.base}"
        OutDict[key]=[SNV(chrom, POS, REF, ALT, line)]

    with open(Calls, 'r') as callsFile:
        linesCalls = [l for l in callsFile if not l.startswith('#') and not l.startswith('"##')]
    for line in linesCalls:
        split=line.split()
        chrom=split[0].split('r')[1]
        if chrom=='X':
            chrom=23
        elif chrom=='Y':
            chrom=24
        else:
            chrom=int(chrom)
        POS=int(split[1])
        REF=split[3]
        ALT=split[4].split(',')
        if "<*>" in ALT:
            ALT.remove("<*>")
        if len(ALT)==0:
            print(line)
            continue
        if split[7].split(';')[0]=="INDEL":
            continue
        depth=int(split[9].split(':')[1])
        AFs=split[9].split(':')[2].split(',')
        for i in range(len(ALT)):
            key=f"{chrom}_{POS}_{REF}_{ALT[i]}"
            currenltAllele=Allele(ALT[i], int(AFs[i+1]), depth)
            if key in OutDict:
                OutDict[key].append(SNV(chrom, POS, REF, currenltAllele, line))

    return(OutDict)

def AFChecker(SNVs, lower, upper, VAF):
    obs=[]
    for key in SNVs:

        truthvar=SNVs[key][0]
        if VAF != None:
            if truthvar.ALT.count != VAF:
                continue

        if len(SNVs[key])==1:
            obs.append(0)
            continue
        
        pilevar=SNVs[key][1]
        obs.append(pilevar.ALT.count/pilevar.ALT.depth)
    print(f"{len(obs)} variants observed")
    print(np.percentile(obs, [lower, upper]))

def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--truthSet', required=True, help='truthset SNVs')
    parser.add_argument('-p', '--pileup', required=True, help='mPileUp Output')
    parser.add_argument('-v', '--vaf', type=float, required=False, help='vaf of variants to look at')
    parser.add_argument('--percentile', type=int, required=False, default=95, help='CI size')
    args = parser.parse_args()
    print("making SNVs")
    SNVs=makeSNVs(args.truthSet, args.pileup)
    lower=(100-args.percentile)/2
    upper=100-lower
    print("comparing Results")
    AFChecker(SNVs, lower, upper, args.vaf)

if __name__ == '__main__':
    main()
