# Indel and SV Alignment Validation

This repository goes over the methodology for validating the Indels and SVs within our somatic variant truthset

## File Structure
```markdown
└── code/
    ├── CigarCaller.py
    ├── CigarCaller.sh
    ├── CigarCallerEX.sh
    └── CombineCigarCalls.sh
```

## Software requirement
* python 3.8 or later with the following packages
  *   argparse
  *   os
  *   collections
* bcftools v1.10.2 or later
* bgzip
* tabix
* truvari

## Overall Pipeline
This pipeline takes in a BAM file and a list of known Indels or SVs as inputs. It will then, for each variant listed, subset the BAM file into just reads overlapping that variant and store them in SAM files. For each read in each of the SAM files, the CIGAR codes will be evaluated to find all insertions or deletions present in any of the reads. Using these CIGAR codes we make a VCF of all the variants present in any of these SAM files, giving a comprehesive list of insertions and deletions. This list is then used as an input to `truvari` to validated the inital list of indels and SVs.

### CigarCaller.sh
The entry point to this pipeline sets up all the directories and files needed for the run. To save time and resources we check the variants on each chromosome seperatly, so this file is used to loop over all the chromosomes

To run this code :

```bash
bash CigarCaller.sh $variant_type $BAM_file $reference_file $sequecing_type $run_name $list_of_variants $out_directory
```
where 
* `variant_type`: "Indel" if making of list of indels or "SV" if making a list of SVs
* `BAM_file`: a BAM file to search for variants in
* `reference_file`: the reference file the BAM file is aligned to
* `sequecing_type`: "True" if the BAM file is long read sequencing, "False" if it is short read sequencing
* `run_name`: a name for the run (will make a directory in the output directory with this name to store the results
* `list_of_variants`: A list of indels or svs to run the pileup on
* `out_directory`: the output directroy to store the results

This code then iterate though each chromosome, and calls ```CigarCallerEX.sh```

### CigarCallerEX.sh
This code subsets the inputed BAM file into various SAM files, one for all the reads overlapping each variant. After creating the SAM files it calls `CigarCaller.py` to call all the indels or SVs in each of the SAM files.

### CigarCaller.py
This python script iterates though each SAM file provided. For each one it will check the CIGAR code for each read, looking for either 'I' or 'D'. When it finds one, it uses the reference fasta along with the read sequence to create the indel/SV and stores it in a dictionary. After iterating though each SAM file, the dictionary represents all variants present. Each of those variants is then written to the output VCF.

### CombineCigarCalls.sh
After the CIGAR calling has been run there will be 24 different VCF outputs, one for each chromosome. The following script combines them all into a single VCF, as well as sorts and indexes the VCF.

To run this code,
```bash
bash CombineCigarCalls.sh $CigarCalls_directory_output $output_vcf_name
```
where ```$CigarCalls_directory_output``` is ```$out_directory\$run_name``` and `output_vcf_name` is the name for the output VCF.

This will produce the file ```$output_directory\$run_name\${output_vcf_name}_sorted.vcf``` as the final VCF output for the entire pipeline.

## Validating the truthset

To use ```truvari``` to validate the truthset, run the following steps:

1. Make sure you are in the correct directory
```bash
cd $out_directory\$run_name
```
2. Run truvari usinging the following paramters
* `-s 1` if indels or `50` if SVs
* `-S 1` if indels or `50` if SVs
* `--pick multi`
* `--pctseq 1` if indels or `0.7` if SVs
* `--pctsize 1` if indels or `0.7` if SVs

Exampe run for SVs:
```bash
truvari bench -b truthset.vcf.gz -c cigarcalls.vcf.gz -f referenceFasta.fna --bSAMple HapMap_Mixture --cSAMple syndip -s 50 -S 50 --pick multi --pctseq 0.7 --pctsize 0.7 --includebed benchmarkRegions.bed -o SV_Validation_Truvari
```
