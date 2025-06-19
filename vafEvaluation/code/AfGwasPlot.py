import argparse
import numpy as np
from scipy.stats import binom
from scipy.stats import binomtest
import sys
import os
import pandas as pd
import math
import matplotlib.pyplot as plt

class Allele:
    def __init__ (self, base, count, depth):
        self.base=base
        self.count=count
        self.depth=depth

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
        return f"chr{self.chrom}:{self.POS}, {self.REF}->{self.ALT.base}"
    
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
        AFString=split[9].split(':')
        if len(AFString)==1:
            AF==0
        elif AFString[1]=='.':
            AF==0
        else:
            AF=float(AFString[1])
        ALT=Allele(ALTBase, AF, 1)
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
            Alleles.append(Allele(ALT[i], int(AFs[i+1]), depth))
        Out.append(SNV(chrom, POS, REF, Alleles, line))
    Out.sort()
    return(Out)

def ReadDepthFile(depthFile):
    dic = {}
    with open(depthFile, "r") as file:
        for line in file:
            split=line.strip().split()
            index=f"{split[0]}_{split[1]}"
            depth=int(split[2])
            dic[index]=depth
    return(dic)

def writeToOutFile(File, Lines):
        f=open(File, "w")
        for line in Lines:
            f.write(f"{line.line}")
        f.close()

def updateStat(n, p, observed):
    if p==0:
        print("error in p value: its 0")
    if n==0:
        observedprob=0
    else:
        observedprob=observed/n
    observedStat=binom.pmf(observed, n, observedprob)
    expectedStat=binom.pmf(observed, n, p)
    if observedStat==0:
        observedStat=sys.float_info.min
    if expectedStat==0:
        expectedStat=sys.float_info.min
    observedTotal=math.log(observedStat)
    expectedTotal=math.log(expectedStat)
    return(observedTotal, expectedTotal)

def makeNullDist(var, depth, iterations):
    Null=[]
    for i in range(iterations):
        p=var.ALT.count
        key=f"{var.chrom}_{var.POS}"
        n=depth
        sample=np.random.binomial(n, p, 1)
        observedTotal, expectedTotal=updateStat(n, p, sample[0])
        if observedTotal==expectedTotal and observedTotal==0:
            return(None)
        Null.append(observedTotal/expectedTotal)
    return Null

def makeChromGWAS(gwas_data, outfile, cuttoff):
    gwas_data['-log10P'] = -np.log10(gwas_data['Pvalue1'])
    gwas_data = gwas_data.sort_values(['CHR', 'POS'])

    plt.figure(figsize=(12, 6))
    colors = ['blue', 'orange']  # Alternating colors for chromosomes

    for i, (chr, group) in enumerate(gwas_data.groupby('CHR')):
        plt.scatter(group['POS'], group['-log10P'], 
                    alpha=0.2, c=colors[i % len(colors)], s=2, label=f'Chromosome {chr}')

    significance_threshold = -np.log10(cuttoff)  # Common GWAS threshold
    plt.axhline(y=significance_threshold, color='red', linestyle='--', label='Significance threshold')
    plt.xlabel("Chromosome")
    plt.ylabel("-log10 p-value")
    plt.title("Manhattan Plot")
    #plt.xticks(chr_boundaries - chr_boundaries.diff().fillna(chr_boundaries.iloc[0]) / 2, range(1, len(chr_boundaries) + 1))
    tickcount=int(len(gwas_data)/20)
    plt.xticks(gwas_data['POS'][::tickcount])
    #plt.legend()
    plt.tight_layout()

    plt.savefig(outfile, dpi=300, format='svg')


def AFChecker(truthSet, mPileup):
    i=0
    counter_1=0
    total=0
    chromcount=0
    cuttoff=0.05/len(truthSet)
    columnNames=["CHR", "POS", "Pvalue1"]
    gwas_overall_data = pd.DataFrame(columns=columnNames)
    gwas_chrom_data = pd.DataFrame(columns=columnNames)
    currentchrom="chr1"
    print(currentchrom)
    for i in range (len(truthSet)):
        truthvar=truthSet[i]
        pilevar=mPileup[i]
        if truthvar.ALT.count==0:
            continue
        if currentchrom!=f"chr{truthvar.chrom}":
            makeChromGWAS(gwas_chrom_data, f"/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SNV/AFPlots/GWAS/AFGWAS{currentchrom}.svg", 0.05/chromcount)
            currentchrom=f"chr{truthvar.chrom}"
            print(currentchrom)
            chromcount=0
            gwas_chrom_data = pd.DataFrame(columns=columnNames)
        p=truthvar.ALT.count
        total+=1
        chromcount+=1
        for alt in pilevar.ALT:
            if alt.base==truthvar.ALT.base:
                n=alt.depth
                p_value = binomtest(alt.count, n, p, alternative='two-sided').pvalue
                if p_value<=cuttoff:
                    counter_1+=1
                gwas_overall_data.loc[len(gwas_overall_data)] = [truthvar.chrom, truthvar.POS, p_value]
                gwas_chrom_data.loc[len(gwas_chrom_data)] = [truthvar.chrom, truthvar.POS, p_value]

    print(f"total variants: {total}")
    print(f"validated by binom method: {total-counter_1} or {(total-counter_1)/total}")

    #gwas_overall_data['-log10P'] = -np.log10(gwas_overall_data['Pvalue1'])
    #gwas_overall_data = gwas_overall_data.sort_values(['CHR', 'POS'])
    #gwas_overall_data['pos_cumulative'] = gwas_overall_data.groupby('CHR')['POS'].cumsum()
    #chr_boundaries = gwas_overall_data.groupby('CHR')['pos_cumulative'].max().cumsum()

    #plt.figure(figsize=(12, 6))
    #colors = ['blue', 'orange']  # Alternating colors for chromosomes

    #for i, (chr, group) in enumerate(gwas_overall_data.groupby('CHR')):
    #    plt.scatter(group['pos_cumulative'], group['-log10P'], 
    #                alpha=0.2, c=colors[i % len(colors)], s=2, label=f'Chromosome {chr}')

    #significance_threshold = -np.log10(cuttoff)  # Common GWAS threshold
    #plt.axhline(y=significance_threshold, color='red', linestyle='--', label='Significance threshold')
    #plt.xlabel("Chromosome")
    #plt.ylabel("-log10 p-value")
    #plt.title("Manhattan Plot")
    #plt.xticks(chr_boundaries - chr_boundaries.diff().fillna(chr_boundaries.iloc[0]) / 2, range(1, len(chr_boundaries) + 1))
    #plt.legend()
    #plt.tight_layout()

    #plt.savefig("/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SNV/AlleleFreqTest/GWASPlotsJan1/AFGWASAll.svg", dpi=300, format='svg')

    #plt.clf()

def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--truthSet', required=True, help='truthset SNVs')
    parser.add_argument('-p', '--pileup', required=True, help='mPileUp Output')
    args = parser.parse_args()
    print("making truth set")
    TruthSetSNVS=makeSNVsFromTruthset(args.truthSet)
    print("making mPileUp set")
    PileUpSNVS=makeSNVsFromPileup(args.pileup)
    print("comparing Results")
    AFChecker(TruthSetSNVS, PileUpSNVS)

    


if __name__ == '__main__':
    main()




        