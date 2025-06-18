# Assembly based Validation

This repository goes over the methodology for validating SNVs, Indels, and SVs usingm an assembly based approach 

## File Structure
```markdown
└── code/
    ├── RegeionReconstruction.py
    ├── RegeionReconstruction.sh
```

## Software requirement
* python 3.8 or later with the following packages
  *   argparse
  *   os
  *   pysam
  *   re
* bcftools v1.10.2 or later
* bgzip
* tabix
* truvari

## Overall Pipeline
This pipeline takes in a list of false positive variants as an input, as well as the path to a vcf for the graph and a path to the directory of the dipcall files (one for each haplotype). The code then works in three steps. First it takes the list of false negitives and splits them by haplotype. So for each haplotype there will be a vcf of all the false negitves from that cell line. In step two the code reconstructs the sequence around each false negitive from each haplotype. It does this by finding all the regions where both the graph and dipcall agree there are no variants (a clean region). Then for each false negitive it find the closest clean region before and after it, and from those regions reconstucts the sequence between twice, first according to the graph, second acourding to dipcall. If these reconstructed sequences as the same then all variants inbetween are validated and outputed to a VCF. For the last step the VCFs for each haplotype are combined into a single VCF

### RegeionReconstruction.sh
The entry point to this pipeline sets up all the directories and files needed for the run. 

to run this code 
```bash
bash CigarCaller.sh $false_neg_vcf $output_dir
```
where 
* false_neg_vcf: the list of false negitive variants
* output_dir: the output directroy to store the results

This code then iterate though each chromosome, and for each one calls ```CigarCallerEX.sh```

### RegeionReconstruction.py
the python script that runs the region reconstuction. It takes in the false negtive set as an input, as well as the path to the graph, dipcall vcfs, and reference file. There is only one paramter to control ``--distance`` which is the minimum size for a clean region. The default value, which was used for this paper, is 200bp


