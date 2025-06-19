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
    
def makeSNVs(truthset, Calls):
    OutDict={}

    with open(truthset, 'r') as truthfile:
        linesTruth = [l for l in truthfile if not l.startswith('#') and not l.startswith('"##')]
    for line in linesTruth:
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
        key=f"{chrom}_{POS}"
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
            #print(line)
            continue
        #print(split[9])
        depth=int(split[9].split(':')[1])
        #print(depth)
        AFs=split[9].split(':')[2].split(',')
        Alleles=[]
        for i in range(len(ALT)):
            Alleles.append(Allele(ALT[i], int(AFs[i+1]), depth))
        key=f"{chrom}_{POS}"
        if key in OutDict:
            OutDict[key].append(SNV(chrom, POS, REF, Alleles, line))

    return(OutDict)

def makeHist(arr, name):
    median = np.median(arr)
    print(median)
    ticks = np.arange(0, 1.1, 0.1)
    bins = np.arange(0, 1.01, 0.01)
    height=0.6087
    width=3.8944

    fig, ax = plt.subplots(figsize=(width, height)) 
    ax.hist(arr, bins=bins, edgecolor='black')
    ax.set_xlabel('Observed Allele Frequnecy', fontdict={'family': 'nimbus sans', 'size': 5}, labelpad=1)
    ax.set_ylabel('Frequency', fontdict={'family': 'nimbus sans', 'size': 5}, labelpad=1)
    ax.set_title(name, fontdict={'family': 'nimbus sans', 'size': 5}, pad=2)
    ax.set_xticks(ticks)
    ax.axvline(median, color='red', linestyle='dashed', linewidth=1, label=f'Median: {median:.4f}')
    ymin, ymax = ax.get_ylim()
    halfway_y = (ymax + ymin) / 2
    ax.text(median + 0.05, halfway_y, f'Median: {median:.4f}', color='black', fontdict={'family': 'nimbus sans', 'size': 5})
    ax.tick_params(axis='x', labelsize=5) 
    ax.tick_params(axis='y', labelsize=5)

    # Save the histogram to a file
    plt.subplots_adjust(top=0.85, bottom=0.40)
    bottom=0.15
    plt.savefig(f"/storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SNV/AFPlots/Hist/CellLineObserved/{name}.svg", format='svg')
    plt.close()
    print(f"graph at /storage1/fs1/jin810/Active/testing/Ruttenberg/SMAHT/SNV/AFPlots/Hist/CellLineObserved/{name}.svg")

def AFChecker(SNVs, outfile):
    obs=[]
    for key in SNVs:
        if len(SNVs[key])==1:
            obs.append(0)
            continue
        
        truthvar=SNVs[key][0]
        pilevar=SNVs[key][1]
        
        for alt in pilevar.ALT:
            if alt.base==truthvar.ALT.base:
                obs.append(alt.count/alt.depth)
    makeHist(obs, outfile)

def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--truthSet', required=True, help='truthset SNVs')
    parser.add_argument('-p', '--pileup', required=True, help='mPileUp Output')
    parser.add_argument('-o', '--Out', required=True, help='outputfile')
    args = parser.parse_args()
    print("making SNVs")
    SNVs=makeSNVs(args.truthSet, args.pileup)
    print(len(SNVs))
    print("comparing Results")
    AFChecker(SNVs, args.Out)


if __name__ == '__main__':
    main()
