# Short Read simulation pipeline

This repository goes over the methodology for simulating the short read sequencing

## File Structure
```markdown
└── cromwell/wdl
    ├── Alignment.wdl
    ├── AlignmentTasks.wdl
    ├── BQSR.wdl
    ├── BwaAlignment.wdl
    ├── LocalRealignment.wdl
    ├── MarkDuplicates.wdl
    ├── RemovePolyGArtifacts.wdl
    ├── RunHapMap.wdl
    ├── Structs.wdl
```

## Running Cromwell workflows
This repositionry will focus on the methodology of simulating short read sequncing and thus will not go over how to set up and run cromwell, as that is system specfic. For more infomation for running cromwells go to LINK.

## Overall Pipeline
This pipeline takes the 10 somatic haplotype fastas (HG002, HG00438, HG02257, HG02486, and HG02622) and simulated short read sequnecing at a specfied depth. Then, for haplotype, the simulated reads get aligned to a reference file. Finailly, these 5 cell lines get combined into a single bam file, which gets merged with a premade HG005 bam file to create the final outputted simulated data.

## RunHapMap.wdl
The entry point to the cromwell pipeline. This script takes the 5 cell lines and generates the maternal and paternal simulation. Then, for each simulation it call ```Alignment.wdl``` to aligned the simulated reads to the reference genome. After all cell lines have been simlated and aligned, this code also merges them together, as well as merges them with HG005. 


Inputs:
* HG002 maternal and paternal fasta
* HG00438 maternal and paternal fasta
* HG02257 maternal and paternal fasta
* HG02486 maternal and paternal fasta
* HG02622 maternal and paternal fasta
* HG005 presimulated bam file 
  
Output:
* HapMap simulated bam file
* HapMap simulated bai file

### wgsim
As a part of this step wgsim is run to generate the simulated reads. It is run with the following parameters:
* ```-N```: the number of reads across the genome to simulated. This was determed based on the read length and read depth desired, as well as the length of the reference genome
* ```-1```: the read length of the first simualted read. For this project 151 was used.
* ```-2```: the read length of the second simualted read. For this project 151 was used.
* ```-r```: the mutation rate (rate of mutation being added to the fasta). For this project we set this to 0.
* ```-e```: the error rate (the rate of a base being misread). For this project we set it to 0 but had it be a controlable paramater.
* ```-R```: the fraction of indels. We set this to 0
* ```-X```: the probability an indel is extended. We set this to 0
* ```-d```: the mean distance between read1 and read2. We set this to 345 to match our real date
* ```-s```: the standard deviation for distance between read1 and read2. We set this to 55 to match our real date


## Alignment.wdl
This wdl script take in a fasta file of simulated reads, and aligns them to a given reference file. This happens in five steps:
* remove polyG artifacts
* alignment with BWA
* mark duplicates
* local realignment
* BQSR

### RemovePolyGArtifacts.wdl
Removes polyG artifact from reads. This is a common artifact from real sequncing data, but not simulated data. However to match the real pipeline for aligning reads, we will perfrom this step in case it intoduces any bias that would be seen in real data

Inputs:
* fasta of simulated reads
  
Output:
* fasta of simulated reads without polyG artifacts

### BwaAlignment.wdl
Use BWA to algin a given set of simulated reads to a refeence genome

Inputs:
* fasta of simulated reads
* refence genome
* index information for reference genome (genetated using BWA index)
  
Output:
* bam file of reads aligned to the reference genome

### MarkDuplicates.wdl
Identified duplicated reads and removes them from the simulation

Inputs:
* bam file of aligned simulated reads
  
Output:
* bam file of aligned simulated reads without duplicated reads

### LocalRealignment.wdl
Indeifed indels in repeat regions and realignes them to standardize the way they are represented (left aligns them)

Inputs:
* bam file of aligned simulated reads
* refence genome
* list of known Indels
  
Output:
* bam file of aligned simulated reads with reads realigned
  
### BQSR.wdl
Recalualex base quality score

Inputs:
* bam file of aligned simulated reads
* refence genome
* list of known SNVs and Indels
  
Output:
* bam file of reads aligned to the reference genome

## Generating HG005 Bam
this pipeline takes in HG005 simulated reads as an input so it doesn't need to be regenerated each time. To generate HG005 we ran the following code

```bash
genomeLenght=3099922541
ReadDepth=417.5
numReads=`awk -v LOG=$genomeLenght 'BEGIN{print (~{ReadDepth}*LOG)/(~{ReadLength}*2*2)}' | awk '{ printf "%.0f\n", $1 }'`
wgsim -N$numReads -1151 -2151} -r0 -e0 -R0 -X0 -d345 -s55 ~{reference} HG005_read1.fa HG005_read2.fa
```
After running wgsim we ran ```alignment.wdl``` to aligned it to HG38
  
