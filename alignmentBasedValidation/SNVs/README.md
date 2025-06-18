# SNV Alignment Validation

This repository goes over the methodology for validating the SNVs within our somatic variant truthset

## File Structure
```markdown
└── code/
    ├── combineMpileup.sh
    ├── runMpileup.sh
    ├── runMpileupEX.sh
└── references
    ├──readgroup.txt

```

## Software requirement
* bcftools v1.10.2 or later
* bgzip
* tabix

## Overall Pipeline
The basis of the pipeline is to run ```bcftools mPileup``` on a bam file to get a list of all variants present within the bam file. We then use ```bcftools isec``` to compare the pileup result with the truthset, giving us our validation rate

### runMpileup.sh
This code sets up all the directories and files needed for the mpileup run. To save time and resources we run mpilelup on each chromosome seperatly, so this file is used to loop over all the chromosomes

to run this code 
```bash
bash runMpileup.sh $bam_file $ref_file $run_name $output_directory $variants_of_interest"
```
where 
* bam_file: the bam file to run the pileup on
* ref_file: the reference file the bam file is aligned to
* run_name: a name for the run (will make a directory in the output directory with this name to store the results
* output_directory: the ooutput directroy to store the results
* variants_of_interest: A list of SNVs to run the pileup on. This is an optional input, if not provided will run a pileup on the entire genome

This code then iterate though each chromosome, and for each one calls ```runMpileupEX.sh```

### runMpileup.sh
this code runs the mpile for the given chromosome. First it takes the inputed variants_of_interest (if provided) and subsets it to the given chromosome. Then it runs mpileup on that chromosome

the following paramaters are used for the mpileup run
* -a FORMAT/AD,FORMAT/DP: formats the output so we can get a count of how many reads show each alternative allele
* --no-BAQ: prevents recomputation of the base alignment quality. This minimizes the number of reads that gets filtered out
* -d 20000: checks the first 20000 reads at a given position, then stops. This makes sure no reads are missed at the position
* -q 0: no map quality filtering
* -Q 0: no base quality filtering
* -G readgroup.txt: tells mpileup that even though the sample is an artifical sample made of combined bam files, treat it as one single bam file
* -R chr#: only run mpileup for the given chromosome

### combineMpileup.sh
After mpileup has been run there will be 24 different vcf outputs, one for each chromosome. The following script combines them all into a single vcf, as well as sorts and indexes the vcf

to run this code 
```bash
bash combineMpileup.sh $mpileup_directory_output
```
where ```$mpileup_directory_output``` is ```$output_directory\$run_name```

this will produce the file ```$output_directory\$run_name\${run_name}_sorted.vcf``` as the final vcf for the entire mpileup pipeline

## validating the truthset

to run ```bcftools isec``` to validate the truthset run the following steps

1. make sure you are in the correct directory
```bash
cd $output_directory\$run_name
```
2. run a multiallelic split on the output vcf. This makes sure when there are multiple variants in the same position they are both matched
```bash
bcftools norm -m - ${run_name}_sorted.vcf -Ov -o ${run_name}_sorted.split.vcf
```
3. use ```bcftools isec``` to compare the variants in the bam file to the truthset
```bash
bcftools isec -p SnvTruthsetComparison $truthset_vcf ${run_name}_sorted.split.vcf
```
This will create a directory ```SnvTruthsetComparison``` to store the comparison.\
The recall will then be the number of variants in ```SnvTruthsetComparison/0002.vcf``` divided by number of variants in the truthset

