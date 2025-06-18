#Docker: biocontainers/python3-pysam:v0.10.0ds-2-deb_cv1
import argparse
import pysam
import os
import re

class variant:
    def __init__ (self, chrom, start, stop, ref, alt):
        self.chrom=chrom
        self.start=start
        self.stop=stop
        self.ref=ref
        self.alt=alt
    def __str__ (self):
        return self.chrom + " " + str(self.start) + "-" + str(self.stop) + " " + self.ref + "->" + self.alt

def makeFastaDict(fasta):
    """
    Input:
        fasta: a reference fasta.
        
    Output:
        FastaDict: A dictoranry where each key is the contig and the entry is the sequences of the contig
    """
    with open(fasta, 'r') as fa:
        lines = [l for l in fa]
    firstLine=lines.pop(0)
    chrom=firstLine.split()[0][1:]
    sequence=""
    FastaDict = {}
    for line in lines:
        if line.startswith(">"):
            FastaDict[chrom]=sequence
            chrom=line.split()[0][1:]
            sequence=""
        else:
            sequence=sequence+line.replace('\n', '')
    FastaDict[chrom]=sequence
    return(FastaDict)

def subsetGraphVariants(GraphVarsList, Sample):
    """
    Input:
        GraphVarsList: a list of variants in the graph
        Sample: the sample the variants are from
        
    Output:
        out: the list of variants convered to a different, more workable form
    """
    out=[]
    for record in GraphVarsList:
        GTIndex=record.samples[Sample]['GT'][0]
        if GTIndex==0 or GTIndex==None:
            continue
        else:
            index=GTIndex-1
        alt=record.alts[index]
        out.append(variant(record.chrom, record.start+1, record.stop, record.ref, alt))
    return(out)

def subsetDipVariants(DipVarList, Sample):
    """
    Input:
        DipVarList: a list of variants in dipcall
        Sample: the sample the variants are from
        
    Output:
        out: the list of variants convered to a different, more workable form
    """
    out=[]
    for record in DipVarList:
        GTIndex=record.samples[Sample]['GT'][1]
        if GTIndex==0 or GTIndex==None:
            continue
        else:
            index=GTIndex-1
        alt=record.alts[index]
        if alt=="*":
            continue
        out.append(variant(record.chrom, record.start+1, record.stop, record.ref, alt))
    return(out)


def splitVars(Negitve, Graph, Dip, reference, graphSample, dipSample, outpath, dist):
    """
    Input:
        Negitve: a list of false negitive variants.
        Graph: the list of variants according to the pangenome graph
        Dip: the list of variants accorsing to dipcall
        reference: the reference fasta for the false negitive variants
        graphSample: the name of the sample in the graph vcf
        dipSample: the name of the sample in the dipcall vcf
        outpath: output directory
        dist: minimum distance between variants in the graph/dipcall. Used to idetfy breakpoints for region reconstruction
    """
    truthList=[]
    FalseList=[]
    region=[]
    MissedVCF = pysam.VariantFile(Negitve)
    GraphVCF = pysam.VariantFile(Graph)
    DipVCF = pysam.VariantFile(Dip)
    header = MissedVCF.header
    for chrom in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]:#, "chrX", "chrY"]:
        print(chrom)
        GraphVar=subsetGraphVariants(list(GraphVCF.fetch(chrom)), graphSample)
        if len(GraphVar)==0:
            continue
        graphFirst=GraphVar[0]
        graphLast=GraphVar[-1]
        GraphCleanRegions=[]
        for i in range(len(GraphVar)-1):
            distance=GraphVar[i+1].start-GraphVar[i].stop-1
            if distance>=dist:
                GraphCleanRegions.append([GraphVar[i].stop+1, GraphVar[i+1].start-1])
        DipVar=subsetDipVariants(list(DipVCF.fetch(chrom)), dipSample)
        DipFirst=DipVar[0]
        DipLast=DipVar[-1]
        DipCleanRegions=[]
        for i in range(len(DipVar)-1):
            distance=DipVar[i+1].start-DipVar[i].stop-1
            if distance>=dist:
                DipCleanRegions.append([DipVar[i].stop+1, DipVar[i+1].start-1])     
        cleanBoth=find_overlaps_linear(GraphCleanRegions, DipCleanRegions, dist)
        firstvar=min(graphFirst.start, DipFirst.start)-1
        lastvar=max(graphLast.stop, DipLast.stop)+1
        cleanRegions=[[firstvar-dist, firstvar]] + cleanBoth + [[lastvar, lastvar+dist]]

        negitiveVars=list(MissedVCF.fetch(chrom))
        i=0
        j=0
        while i<len(negitiveVars):
            var=negitiveVars[i]
            #print(f"start it: {var}")
            for record in DipVCF.fetch(var.chrom, var.pos-1, var.pos):  # Position is 1-based, VCF is 0-based
            # Check if the position, reference, and alternate allele match
                if record.chrom == var.chrom and record.pos == var.pos and record.ref == var.ref and var.alts[0] in record.alts[0]:
                    #print(f"exact match for {var}")
                    truthList.append(var)
                    i+=1
                    break
            beforeRegion=cleanRegions[j]
            AfterRegion=cleanRegions[j+1]
            if beforeRegion[1] >= var.start+1 or AfterRegion[0] <= var.stop:
                j+=1
            else:
                start=beforeRegion[1]-100
                end=AfterRegion[0]+100
                if compareRegions(reference[chrom], chrom, start, end, GraphVCF, DipVCF, graphSample, dipSample):
                    truthList.append(var)
                else:
                    FalseList.append(var)
                i+=1
    #print("validated Variants")
    #print(len(truthList))
    #print("unvalidated Variants")
    #print(len(FalseList))
    FalsePath=outpath + "/fn.vcf"
    TruePath=outpath + "/tp.vcf"
    f=open(FalsePath, "w")
    f.write(str(header))
    for record in FalseList:
        f.write(str(record))
    f.close()
    t=open(TruePath, "w")
    t.write(str(header))
    for record in truthList:
        t.write(str(record))
    t.close()

def find_overlaps_linear(list1, list2, dist):
    """
    Input:
        list1: a list of genomic regions
        list2: a second list of genomic regions
        dist: minimum distance for an output region to be kept
    Output:
        overlaps: a list of regions in both list1 and list2
    """
    overlaps = []
    i, j = 0, 0 

    while i < len(list1) and j < len(list2):
        start1, end1 = list1[i]
        start2, end2 = list2[j]

        # Calculate the overlap
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)

        # If there's a valid overlap, add it to the result
        if overlap_start <= overlap_end and overlap_end-overlap_start>=dist:
            overlaps.append([overlap_start, overlap_end])

        # Move the pointer for the region that ends first
        if end1 < end2:
            i += 1
        else:
            j += 1

    return overlaps


def makeAssembly(Variants, Reference, chrom, start, end):
    """
    Input:
        Variants: a list of variants in a region
        Sample: the sample the variants are from
        chrom: the chromosome of that region
        start: the start position of the region
        end: the end position of the region
        
    Output:
        seq: the recontsuctred sequence of the region
    """
    pos=start
    seq=""
    for var in Variants:
        seq+=Reference[pos-1:var.start-1]
        seq+=var.alt
        pos=var.start+len(var.ref)
    seq+=Reference[pos-1:end-1]
    return(seq)


def compareRegions(Reference, chromosome, start, end, GraphVCF, DipcallVCF, GraphSample, DipSample):
    """
    Input:
        Reference: the reference file
        chrom: the chromosome of a region
        start: the start position of the region
        end: the end position of the region
        GraphVCF: the list of variants according to the pangenome graph
        DipcallVCF: the list of variants accorsing to dipcall
        GraphSample: the name of the sample in the graph vcf
        DipSample: the name of the sample in the dipcall vcf
    Output:
        a boolean on if the sequnece of the region is the same in dipcall compared to the graph
    """
    GraphVars=subsetGraphVariants(GraphVCF.fetch(chromosome, start, end), GraphSample)
    DipVars=subsetDipVariants(DipcallVCF.fetch(chromosome, start, end), DipSample)
    GraphRegion=makeAssembly(GraphVars, Reference, chromosome, start, end)
    DipRegion=makeAssembly(DipVars, Reference, chromosome, start, end)
    #if (GraphRegion!=DipRegion):
    #    print(f"mismatch: {chromosome}:{start}-{end}")
    #    print(f"graph: {GraphRegion}")
    #    print(f"dip: {DipRegion}")
    return(GraphRegion==DipRegion)

def main(): 
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--Missed', required=True, help='VCF Missed Variants to rescue')
    parser.add_argument('-g', '--Graph', required=True, help='VCF of Graph Variants')
    parser.add_argument('-d', '--Dip', required=True, help='VCF of Dipcall Variants')
    parser.add_argument('-f', '--Ref', required=True, help='Reference Fasta')
    parser.add_argument('--GraphSample', required=True, help='Graph Sample Column')
    parser.add_argument('--DipSample', required=True, help='Dipcall Sample Column')
    parser.add_argument('-o', '--outDir', required=True, help='Output directory path')
    parser.add_argument('--distance', type=int, required=False, default=200, help='distance between two graph variants to be considered a reference region')

    args = parser.parse_args()
    os.makedirs(args.outDir, exist_ok=False)
    print("Making reference Dict")
    RefDict=makeFastaDict(args.Ref)
    print("split Vars")
    splitVars(args.Missed, args.Graph, args.Dip, RefDict, args.GraphSample, args.DipSample, args.outDir, args.distance)

if __name__ == '__main__':
    main()



