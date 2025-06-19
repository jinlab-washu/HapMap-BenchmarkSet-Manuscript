# Variant Allele Frqunecy Validation Code

This repository goes over all the code used to validate the VAF of the variants in the truthset

## File Structure
```markdown
└── code/
    ├── AfGwasPlot.py
    ├── expectedVsObservedVafGraph.py
    ├── graphObserevedVAF.py
    ├── percentiles.py
    ├── validateSNVsByVAF.py
```

## AfGwasPlot.py

### Overall Pipeline
This pipeline produces the top image in figure 2e. In summary, the code takes in a vcf of the truthset, and the pileup of those variant.\
Then, for each variant in the truthset, it identies the expected allele frequnecy and observed read depth at the position. Using this is creates\
a binomial model for the observed allele frequnecy and find the probability of the actual observed allele fruency being observed. It then \
plots the value in a large GWAS like plot.

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
python AfGwasPlot.py -t $truthset_vcf -p $mpileup_vcf
```
where ```truthset_vcf``` is a vcf of the truthset variants and ```mpileup_vcf``` is a vcf of the mpileup output for the truthset variants

## expectedVsObservedVafGraph.py

### Overall Pipeline
This pipeline graphs the difference between the expected allele frequency and the observed allele frequency in both a linear and log scale, as well as plotting the expected allele frequnecy vs the observed allele frequency

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * os
  * matplotlib.pyplot
    
### Execution
To execute the code run the following command
```bash
python expectedVsObservedVafGraph.py -t $truthset_vcf -p $mpileup_vcf -o $output_directroy
```
where ```truthset_vcf``` is a vcf of the truthset variants, ```mpileup_vcf``` is a vcf of the mpileup output for the truthset variants, and ```output_directroy``` is the output directory to store the graphs

## graphObserevedVAF.py

### Overall Pipeline
This pipeline graphs the observed allele frequency of a group of variants inputed by the user as a histogram.

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * matplotlib.pyplot

### Execution
To execute the code run the following command
```bash
python graphObserevedVAF.py -t $truthset_vcf -p $mpileup_vcf -o $output_file
```
where ```truthset_vcf``` is a vcf of the truthset variants, ```mpileup_vcf``` is a vcf of the mpileup output for the truthset variants, and ```output_file``` is the path to the file to store the historgram

## percentiles.py

### Overall Pipeline
This pipleine idenfies the lower and upper observed VAF bounds to captrues a given percent of the variants in the truth set. So it 95% is inputed, the code finds the 2.5 percentile VAF and the 97.5 percentile VAF

### Software requirments
* python3.8 with the following packages
  * argparse
  * numpy
  * matplotlib.pyplot

### Execution
To execute the code run the following command
```bash
python percentiles.py -t $truthset_vcf -p $mpileup_vcf -v $VAF --percentile $percentile
```
where 
* ```truthset_vcf```: a vcf of the truthset variants
* ```mpileup_vcf```: a vcf of the mpileup output for the truthset variants
* ```VAF``` an optional float of what VAF to consider. If provided will only look at variants with the given expected VAF
* ```percentile```: what percentile to spand (default: 95)
  
## validateSNVsByVAF.py

### Overall Pipeline
This pipleline looks at all the SNVs validated by alignment based methodolgies and find how many of them have an observed VAF close to the expected VAF, as defined by a binomial distibution. For each variant in the truthset, it identies the expected allele frequnecy and observed read depth at the position. Using this is creates a binomial model for the observed allele frequnecy and find the probability of the actual observed allele fruency being observed. For each variant with a small enough probability (0.05/# of variants) it is considered an unvalidated VAF, otherwise it is validated

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
python validateSNVsByVAF.py -t $truthset_vcf -p $mpileup_vcf
```
where ```truthset_vcf``` is a vcf of the truthset variants and ```mpileup_vcf``` is a vcf of the mpileup output for the truthset variants
