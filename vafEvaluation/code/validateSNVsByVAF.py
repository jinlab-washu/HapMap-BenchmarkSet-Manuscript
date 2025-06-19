import argparse
import numpy as np
from scipy.stats import binomtest
from scipy.stats import chi2_contingency
import sys
import os
import pandas as pd
import math
import re

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
    """
    Input:
        file: the vcf of the truthset variants
        
    Output:
        SNVList: A list of SNVs reformated to be more workable
    """
    SNVList=[]
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
        SNVList.append(SNV(chrom, POS, REF, ALT, line))
    SNVList.sort()
    return(SNVList)

def makeSNVsFromPileup(file):
    """
    Input:
        file: the vcf of the pileup variants
        
    Output:
        SNVList: A list of SNVs reformated to be more workable
    """
    SNVList=[]
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
        SNVList.append(SNV(chrom, POS, REF, Alleles, line))
    SNVList.sort()
    return(SNVList)

def writeToOutFile(File, Lines):
    """
    Input:
        File: path to file to write to
        Lines: Lines to write to the file
    """
    f=open(File, "w")
    for line in Lines:
        f.write(f"{line.line}")
    f.close()

def checkBed (var, bed):
    """
    Input:
        var: a single SNV
        bed: a list of bed files for the disotrion regions in each haplotype
    """
    info=var.line.strip().split()[7]
    genotypes = re.findall(r'[0-9.*]+\|[0-9.*]+', info)
    numbers = [
        int(allele) if allele.isdigit() else 0
        for pair in genotypes
        for allele in pair.split('|')
    ]
    for i in range(0,12):
        if numbers[i]==1:
            with open(bed[i], 'r') as f:
                for line in f:
                    bedsplit=line.strip().split()
                    if var.POS>=int(bedsplit[1]) and var.POS<int(bedsplit[2]):
                        return True
    return False

def AFChecker(truthSet, mPileup, bedList):
    """
    Input:
        truthSet: the truthset SNVs
        mPileup: the pileup SNVs
        bedList: a list of bed files for the disotrion regions in each haplotype
    """    
    i=0
    counter=0
    total=0
    cuttoff=0.05/len(truthSet)
    toWrite=[]

    print(f"cutoff: {cuttoff}")
    currentchrom="chr1"
    print(currentchrom)

    validatedDist=0
    unvalidatedDist=0
    validatedNoDist=0
    unvalidatedNoDist=0

    for i in range (len(truthSet)):
        truthvar=truthSet[i]
        pilevar=mPileup[i]
        if truthvar.ALT.count==0:
            continue
        if currentchrom!=f"chr{truthvar.chrom}":
            currentchrom=f"chr{truthvar.chrom}"
            print(currentchrom)
            break
        p=truthvar.ALT.count
        total+=1
        dist=checkBed(truthvar, bedList)
        for alt in pilevar.ALT:
            if alt.base==truthvar.ALT.base:
                n=alt.depth
                p_value = binomtest(alt.count, n, p, alternative='two-sided').pvalue
                if p_value<=cuttoff:
                    counter+=1
                    if dist:
                        unvalidatedDist+=1
                    else:
                        unvalidatedNoDist+=1
                else:
                    if dist:
                        validatedDist+=1
                    else:
                        validatedNoDist+=1
    data = np.array([[validatedNoDist, validatedDist],
                 [unvalidatedNoDist, unvalidatedDist]])
    chi2, p, dof, expected = chi2_contingency(data)



    print(f"total variants: {total}")
    print(f"validated by binom method: {total-counter} or {(total-counter)/total}")
    print(f"validated distortion: {validatedDist}")
    print(f"unvalidated distortion: {unvalidatedDist}")
    print(f"validated not in distortion: {validatedNoDist}")
    print(f"unvalidated not in distortion: {unvalidatedNoDist}")
    print(f"chi squared p vallue: {p}")

def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--truthSet', required=True, help='truthset SNVs')
    parser.add_argument('-p', '--pileup', required=True, help='mPileUp Output')
    args = parser.parse_args()

    HG002_Pat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG002_Pat_dist.bed"
    HG002_Mat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG002_Mat_dist.bed"
    HG00438_Pat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG00438_Pat_dist.bed"
    HG00438_Mat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG00438_Mat_dist.bed"
    HG005_Pat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG005_Pat_dist.bed"
    HG005_Mat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG005_Mat_dist.bed"
    HG02257_Pat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG02257_Pat_dist.bed"
    HG02257_Mat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG02257_Mat_dist.bed"
    HG02486_Pat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG02486_Pat_dist.bed"
    HG02486_Mat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG02486_Mat_dist.bed"
    HG02622_Pat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG02622_Pat_dist.bed"
    HG02622_Mat_Dist="/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/References/Blackout/Distortion/BlackoutOptions/DistByCellLine/HG02622_Mat_dist.bed"
    bedList=[HG002_Mat_Dist, HG002_Pat_Dist, HG00438_Mat_Dist, HG00438_Pat_Dist, HG005_Mat_Dist, HG005_Pat_Dist, HG02257_Mat_Dist, HG02257_Pat_Dist, HG02486_Mat_Dist, HG02486_Pat_Dist, HG02622_Mat_Dist, HG02622_Pat_Dist]

    print("making truth set")
    TruthSetSNVS=makeSNVsFromTruthset(args.truthSet)
    print("making mPileUp set")
    PileUpSNVS=makeSNVsFromPileup(args.pileup)
    print("comparing Results")
    AFChecker(TruthSetSNVS, PileUpSNVS, bedList)

    


if __name__ == '__main__':
    main()




        
