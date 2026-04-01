# HapMap-BenchmarkSet-Manuscript

This repository holds scripts for all analysis done for [this paper](), including benchmark set validation, model fitting, and figure generation. More information about can be found in Methods section of our paper. 

## File Structure
```markdown
├── alignmentBasedValidation/
    ├── IndelsAndSVs/
    └── SNVs/
├── assemblyBasedValidation/
    └── regionRescue/
├── manuscriptFigures/
├── simulations/
    ├── shortReadSimulation/
    └── longReadSimulation/
└── vafEvaluation/
```


## Short and Long Read Simulation

This repository denotes code used for generating short- and long-read sequencing data that were used in subsequent benchmark set validation steps.

## Alignment Based Validation

This repository goes over the methodology for validating the SNVs, Indels, and SVs within our somatic variant benchmarkset using an alignment based approach. See READMEs in subdirectories for more information.

## Assembly Based Validation

This repository goes over the methodology for validating SNVs, Indels, and SVs using an assembly based approach. See READMEs in subdirectories for more information.

## Variant Allele Frequency Evaluation

This repository goes over all the code used to validate the VAF of the variants in the benchmarkset. See READMEs in subdirectories for more information.

## Manuscript Figure Generation

This repository holds scripts to generate all figures in our manuscript including caller performance analysis, and gene level coverage requirement analysis.


    
