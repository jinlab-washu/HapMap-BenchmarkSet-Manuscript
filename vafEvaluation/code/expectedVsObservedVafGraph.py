import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

class Allele:
    def __init__ (self, base, freq):
        self.base=base
        self.freq=freq

    def __eq__ (self, other):
        return self.base==other.base
    def __str__ (self):
        return f"base: {self.base}, freq: {self.freq}"

class SNV:
    def __init__ (self, chrom, POS, REF, ALT, line):
        self.chrom=chrom
        self.POS=POS
        self.REF=REF
        self.ALT=ALT
        self.line=line
    
    def __str__ (self):
        return f"chr{self.chrom}:{self.POS}"

    def __eq__(self, other):
        return self.chrom==other.chrom and self.POS==other.POS
    
    def __lt__(self, other):
        if self.chrom!=other.chrom:
            return self.chrom<other.chrom
        return self.POS<other.POS

    def compareSNV (self, other):
        if self.chrom!=other.chrom or self.POS!=other.POS:
            return False
        return self.ALT in other.ALT

        
              
    
def makeSNVsFromTruthset(file):
    Out=[]
    with open(file, 'r') as of:
        lines = [l for l in of if not l.startswith('#') and not l.startswith('"##')]
    for line in lines:
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
        ALTBase=split[4]
        AFString=split[9].split(':')[1]
        if AFString=='.':
            AFString='0'
        AF=float(AFString)
        ALT=Allele(ALTBase, AF)
        Out.append(SNV(chrom, POS, REF, ALT, line))
    Out.sort()
    return(Out)

def makeSNVsFromPileup(file):
    Out=[]
    with open(file, 'r') as of:
        lines = [l for l in of if not l.startswith('#') and not l.startswith('"##')]
    for line in lines:
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
            continue
        if split[7].split(';')[0]=="INDEL":
            continue
        depth=int(split[9].split(':')[1])
        AFs=split[9].split(':')[2].split(',')
        Alleles=[]
        for i in range(len(ALT)):
            Alleles.append(Allele(ALT[i], int(AFs[i+1])/depth))
        Out.append(SNV(chrom, POS, REF, Alleles, line))
    Out.sort()
    return(Out)

def writeToOutFile(File, Lines):
        f=open(File, "w")
        for line in Lines:
            if line.chrom==22:
                chromeString='X'
            elif line.chrom==23:
                chromeString='Y'
            else:
                chromeString=line.chrom
            f.write(f"{line.line}")
        f.close()

def makeChart(truthSet, mPileup, outDir):
    expectedAF=[]
    observedAF=[]
    diff=[]
    for i in range(len(truthSet)):
        truthsetVar=truthSet[i]
        pileVar=mPileup[i]
        if truthsetVar!=pileVar:
            print("something just broke")
            break
        if truthsetVar.ALT.freq==0:
            continue
        for alt in pileVar.ALT:
            if alt.base==truthsetVar.ALT.base:
                expectedAF.append(truthsetVar.ALT.freq)
                observedAF.append(alt.freq)
                diff.append(truthsetVar.ALT.freq-alt.freq)
                break

    height=2.0096
    width=2.7431
    sizeAxis=5
    sizeTitle=7

    fig, ax = plt.subplots(figsize=(width, height)) 
    ax.scatter(expectedAF, observedAF, color='#0077BB', marker='o', s=10, alpha=0.3)
    ax.set_xlabel('Expected Allele Frequency', fontdict={'family': 'nimbus sans', 'size': sizeAxis})
    ax.set_ylabel('Observed Allele Frequency', fontdict={'family': 'nimbus sans', 'size': sizeAxis})
    ax.set_title('Expected vs Observed Allele Frequency', fontdict={'family': 'nimbus sans', 'size': sizeTitle})
    ax.tick_params(axis='x', labelsize=sizeAxis) 
    ax.tick_params(axis='y', labelsize=sizeAxis)
    plt.subplots_adjust(bottom=0.20)
    plt.savefig(f"{outDir}/ExpectedVsObservedAF.svg", format='svg')

    plt.clf()

    fig, ax = plt.subplots(figsize=(width, height)) 
    ax.hist(diff, bins=100, color='#0077BB', weights=[100 / len(diff)] * len(diff), alpha=0.7, edgecolor='black')
    ax.set_xlabel('Difference between expected and observed', fontdict={'family': 'nimbus sans', 'size': sizeAxis})
    ax.set_ylabel('Percent', fontdict={'family': 'nimbus sans', 'size': sizeAxis})
    ax.set_title('Expected vs Observed Allele Frequency', fontdict={'family': 'nimbus sans', 'size': sizeTitle})
    ax.tick_params(axis='x', labelsize=sizeAxis) 
    ax.tick_params(axis='y', labelsize=sizeAxis)
    plt.subplots_adjust(bottom=0.20)
    plt.savefig(f"{outDir}/DifferenceInAF.svg", format='svg')

    plt.clf()

    fig, ax = plt.subplots(figsize=(width, height)) 
    ax.hist(diff, bins=100, color='#0077BB', weights=[100 / len(diff)] * len(diff), alpha=0.7, edgecolor='black')
    ax.set_yscale('log')
    ax.set_xlabel('Difference between expected and observed', fontdict={'family': 'nimbus sans', 'size': sizeAxis})
    ax.set_ylabel('Percent (log scale)', fontdict={'family': 'nimbus sans', 'size': sizeAxis})
    ax.set_title('Expected vs Observed Allele Frequency', fontdict={'family': 'nimbus sans', 'size': sizeTitle})
    ax.tick_params(axis='x', labelsize=sizeAxis) 
    ax.tick_params(axis='y', labelsize=sizeAxis)
    plt.subplots_adjust(bottom=0.20)
    plt.savefig(f"{outDir}/DifferenceInAFLog.svg", format='svg')



def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--truthSet', required=True, help='truthset SNVs')
    parser.add_argument('-p', '--pileup', required=True, help='mPileUp Output')
    parser.add_argument('-o', '--outDir', required=True, help='output directory for graphs')
    args = parser.parse_args()
    print("making truth set")
    TruthSetSNVS=makeSNVsFromTruthset(args.truthSet)
    print("making mPileUp set")
    PileUpSNVS=makeSNVsFromPileup(args.pileup)
    print("comparing Results")
    makeChart(TruthSetSNVS, PileUpSNVS, args.outDir)

    


if __name__ == '__main__':
    main()




        
