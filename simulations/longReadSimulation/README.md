# Long-Read Simulation Pipeline

This repository provides a comprehensive methodology for simulating long-read PacBio sequencing data with quality control and analysis workflows.

## File Structure
```markdown
├── example_script/
    ├── step1_pbsim_HG002_mat_split1.bjob
    ├── step2_ccs_HG002_mat_split1_0.bjob
    ├── step2_ccs_HG002_mat_split1_0.sh
    ├── step3_check_ccs_report_HG002.sh
    └── step4_generate_stat_HG002_mat_1_1.bjob
├── ccs_report_check.py
├── main.py
├── bam_stat.py
└── pbsim3p.py
```

## Overview

This pipeline simulates high-coverage PacBio long-read sequencing data for HapMap samples (HG002, HG005, HG00438, HG02257, HG02486, HG02622) at various coverage depths, with comprehensive quality control and statistical analysis.

## Pipeline Workflow

### Main driver script
**Script**: `main.py`

### Step 1: Long-Read Simulation with PBSIM3
- **Script**: `pbsim3p.py`
- **Purpose**: Generate simulated PacBio reads from reference assemblies
- **Key Features**:
  - Parallelized execution using `pbsim3p.py` wrapper
  - Configurable coverage depths
  - Sample-specific coverage ratios based on experimental design
  - Automated job scheduling with LSF

**PBSIM3 Parameters**:
- Strategy: `--strategy wgs` (whole genome sequencing)
- Method: `--method qshmm` with QSHMM model
- Coverage: Variable (40x for HG002, 1670x for HG005, etc.)
- Pass number: `--pass-num 7`
- Read length: Mean 22,000bp, min 300bp, max 80,000bp
- Accuracy: `--accuracy-mean 1`

### Step 2: Consensus Sequence Generation (CCS)
- **Scripts**: `step2_ccs_*.bjob`, `step2_ccs_*.sh`
- **Purpose**: Generate high-quality consensus reads from simulated subreads
- **Process**:
  1. Index BAM files using `pbindex`
  2. Generate CCS reads using PacBio `ccs` tool


### Step 3: Quality Control Validation
- **Scripts**: `ccs_report_check.py`, `step3_check_ccs_report_*.sh`
- **Purpose**: make sure the CCS reads were generated correctly
- **Validation Criteria**: ZMW pass filter rate

### Step 4: Statistical Analysis and Visualization
- **Scripts**: `bam_stat.py`, `step4_generate_stat_*.bjob`
- **Purpose**: Generate read statistics and quality control figures
- **Analyses**:
  - Read length distribution comparison
  - Number of passes distribution
  - GC content analysis
  - Per-read average PHRED score distribution
  - Post-alignment MAPQ distribution
  - Coverage distribution analysis

## Sample Configuration

The pipeline processes six HapMap samples with specified coverage ratios:

| Sample | Coverage Ratio | Target Coverage (per parental haplotype) |
|--------|---------------|------------------------------------------|
| HG005  | 83.5%         | 1670x                                    |
| HG02622| 10.0%         | 200x                                     |
| HG002  | 2.0%          | 40x                                      |
| HG02257| 2.0%          | 40x                                      |
| HG02486| 2.0%          | 40x                                      |
| HG00438| 0.5%          | 10x                                      |

**Total Target Coverage**: 4000x

## Running the Pipeline

### Prerequisites
- PBSIM3 installed and configured
- PacBio tools (pbindex, ccs) available
- Python 3.x with required libraries (matplotlib, seaborn, numpy)

### Execution Steps

1. **Configure Sample Parameters** in `main.py`:
   ```python
   sample_ratios = {
       "HG005": 0.835,
       "HG02622": 0.10,
       # ... other samples
   }
   total_coverage = 4000
   ```

2. **Run PBSIM3 Simulation**:
   ```bash
   python main.py  # Generates step1 job scripts
   # Submit generated .bjob files to LSF
   ```

3. **Execute CCS Generation**:
   ```bash
   # Submit step2 job scripts after step1 completion
   ```

4. **Validate Quality**:
   ```bash
   bash step3_check_ccs_report_*.sh
   ```

5. **Generate Statistics and Figures**:
   ```bash
   # Submit step4 job scripts
   jupyter notebook Supp_Figure1.ipynb
   ```
   Note, for organization purposes, figure generation scripts have been stored in `../manuscriptFigures/Supp_Figure1.ipynb`.

## Output Files

- **Simulated Reads**: BAM files with simulated PacBio subreads
- **CCS Reads**: High-quality consensus sequences
- **Read Summaries**: Read length, quality, and coverage distributions
- **Visualization**: Comparative plots between real and simulated data

## Resource Requirements

### Computational Resources
- **CPU**: 30-31 cores per job (configurable)
- **Memory**: 30-50GB per job depending on step
- **Storage Estimates**
    - 10x coverage: ~166GB compressed
    - 40x coverage: ~655GB compressed  
    - 200x coverage: ~3.2TB compressed

