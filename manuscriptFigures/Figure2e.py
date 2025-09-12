#################################################################
# for every variant, plots the probability of the observed allele
# frequnecy given the null hypothesis that observed allele frequnecy 
# follows a binomial distibution when p=expected allele frqunecy and 
# n=read depth at the specific position

## author: Andrew Ruttenberg
## contact: ruttenberg.andrew@wustl.edu
#################################################################

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
              
    
def makeSNVsFromBenchmarkset(file):
    """
    Input:
        file: the vcf of the benchmarkset variants
        
    Output:
        SNVList: A list of SNVs reformated to be more workable
    """
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
    """
    Input:
        file: the vcf of the pileup variants
        
    Output:
        SNVList: A list of SNVs reformated to be more workable
    """
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
    """
    Input:
        depthFile: a file with the read depth at each position
    Output:
        depthDict: a dictonary where each index is a genomic position and each entry is the read depth
    """
    depthDict = {}
    with open(depthFile, "r") as file:
        for line in file:
            split=line.strip().split()
            index=f"{split[0]}_{split[1]}"
            depth=int(split[2])
            depthDict[index]=depth
    return(depthDict)

def makeChromGWAS(gwas_data, outfile, cuttoff):
    """
    Input:
        gwas_data: a data frame with the read depth, expected VAF, and observed VAF for each variant being plotted
        outfile: path to the output file
        cuttoff: significance threshold 
    """
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


def AFChecker(benchmarkSet, mPileup):
    """
    Input:
        benchmarkSet: list of benchmarkset variants
        outfile: list of pileup variants
    """
    i=0
    counter_1=0
    total=0
    chromcount=0
    cuttoff=0.05/len(benchmarkSet)
    columnNames=["CHR", "POS", "Pvalue1"]
    gwas_overall_data = pd.DataFrame(columns=columnNames)
    gwas_chrom_data = pd.DataFrame(columns=columnNames)
    currentchrom="chr1"
    print(currentchrom)
    for i in range (len(benchmarkSet)):
        benchmarkvar=benchmarkSet[i]
        pilevar=mPileup[i]
        if benchmarkvar.ALT.count==0:
            continue
        if currentchrom!=f"chr{benchmarkvar.chrom}":
            makeChromGWAS(gwas_chrom_data, f"/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SNV/AFPlots/GWAS/AFGWAS{currentchrom}.svg", 0.05/chromcount)
            currentchrom=f"chr{benchmarkvar.chrom}"
            print(currentchrom)
            chromcount=0
            gwas_chrom_data = pd.DataFrame(columns=columnNames)
        p=benchmarkvar.ALT.count
        total+=1
        chromcount+=1
        for alt in pilevar.ALT:
            if alt.base==benchmarkvar.ALT.base:
                n=alt.depth
                p_value = binomtest(alt.count, n, p, alternative='two-sided').pvalue
                if p_value<=cuttoff:
                    counter_1+=1
                gwas_overall_data.loc[len(gwas_overall_data)] = [benchmarkvar.chrom, benchmarkvar.POS, p_value]
                gwas_chrom_data.loc[len(gwas_chrom_data)] = [benchmarkvar.chrom, benchmarkvar.POS, p_value]

    print(f"total variants: {total}")
    print(f"validated by binom method: {total-counter_1} or {(total-counter_1)/total}")

def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--benchmarkSet', required=True, help='benchmarkset SNVs')
    parser.add_argument('-p', '--pileup', required=True, help='mPileUp Output')
    args = parser.parse_args()
    print("making benchmark set")
    BenchmarkSetSNVS=makeSNVsFromBenchmarkset(args.benchmarkSet)
    print("making mPileUp set")
    PileUpSNVS=makeSNVsFromPileup(args.pileup)
    print("comparing Results")
    AFChecker(BenchmarkSetSNVS, PileUpSNVS)

if __name__ == '__main__':
    main()




        
