# Variant Allele Frqunecy Validation Code

This repository goes over all the code used to validate the VAF of the variants in the benchmarkset.

## File Structure
```markdown
└── code/
    ├── AfGwasPlot.py
    ├── expectedVsObservedVafGraph.py
    ├── graphObserevedVAF.py
    ├── percentiles.py
    └── validateSNVsByVAF.py
```

## AfGwasPlot.py

### Summary
This python script produces the top image in figure 2e. In summary, the code takes in a VCF of the benchmarkset, and the pileup of those variant.\
Then, for each variant in the benchmarkset, it identifies the expected allele frequnecy and observed read depth at the position. Using this is creates a binomial model for the observed allele frequnecy and find the probability of the actual observed allele frequency being observed. It then 
plots the values in a large GWAS-like plot.

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * scipy.stats
  * sys
  * os
  * pandas
  * math
  * matplotlib.pyplot   

### Execution
To execute the code run the following command
```bash
python AfGwasPlot.py -t $benchmarkset_vcf -p $mpileup_vcf
```
where ```benchmarkset_vcf``` is a VCF of the benchmarkset variants and ```mpileup_vcf``` is a VCF of the mpileup output for the benchmarkset variants.

## expectedVsObservedVafGraph.py

### Summary
This python script graphs the difference between the expected allele frequency and the observed allele frequency in both linear and log scales, as well as plotting the expected allele frequnecy vs the observed allele frequency.

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * os
  * matplotlib.pyplot
    
### Execution
To execute the code run the following command
```bash
python expectedVsObservedVafGraph.py -t $benchmarkset_vcf -p $mpileup_vcf -o $output_directroy
```
where ```benchmarkset_vcf``` is a VCF of the benchmarkset variants, ```mpileup_vcf``` is a VCF of the mpileup output for the benchmarkset variants, and ```output_directroy``` is the output directory to store the graphs.

## graphObserevedVAF.py

### Summary
This python script graphs the observed allele frequency of a group of variants inputed by the user as a histogram.

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * matplotlib.pyplot

### Execution
To execute the code run the following command
```bash
python graphObserevedVAF.py -t $benchmarkset_vcf -p $mpileup_vcf -o $output_file
```
where ```benchmarkset_vcf``` is a VCF of the benchmarkset variants, ```mpileup_vcf``` is a VCF of the mpileup output for the benchmarkset variants, and ```output_file``` is the path to the file to store the historgram.

## percentiles.py

### Summary
This python script idenfies the lower and upper observed VAF bounds to captrues a given percent of the variants in the benchmark set. So if 95% is inputed, the code finds the 2.5 percentile VAF and the 97.5 percentile VAF.

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * matplotlib.pyplot

### Execution
To execute the code run the following command
```bash
python percentiles.py -t $benchmarkset_vcf -p $mpileup_vcf -v $VAF --percentile $percentile
```
where 
* ```benchmarkset_vcf```: a VCF of the benchmarkset variants
* ```mpileup_vcf```: a VCF of the mpileup output for the benchmarkset variants
* ```VAF``` an optional float of what VAF to consider. If provided will only look at variants with the given expected VAF
* ```percentile```: what percentile to spand (default: 95)
  
## validateSNVsByVAF.py

### Summary
This python script looks at all the SNVs validated by alignment based methodologies and find how many of them have an observed VAF close to the expected VAF, as defined by a binomial distibution. For each variant in the benchmarkset, it identifies the expected allele frequnecy and observed read depth at the position. Using this it creates a binomial model for the observed allele frequnecy and find the probability of the actual observed allele fruency being observed. For each variant with a small enough probability (0.05/# of variants) it is considered an unvalidated VAF, otherwise it is validated.

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * matplotlib.pyplot
  * scipy.stats
  * sys
  * pandas
  * math
  * re
    
### Execution
To execute the code run the following command
```bash
python validateSNVsByVAF.py -t $benchmarkset_vcf -p $mpileup_vcf
```
where ```benchmarkset_vcf``` is a VCF of the benchmarkset variants and ```mpileup_vcf``` is a VCF of the mpileup output for the benchmarkset variants.
