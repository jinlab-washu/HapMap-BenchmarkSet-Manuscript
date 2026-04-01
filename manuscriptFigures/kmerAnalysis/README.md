# kmer-enrichment-analysis

Last updated by Zitian Tang. 2026/04/01.

Identifies sequence composition features that distinguish **benchmark SV-harboring regions** from **regions with no true-positive SV calls (noTP)** in the human genome. The pipeline defines each region set, extracts the corresponding sequences from GRCh38, counts k-mers at k = 31 with Jellyfish, and tests for differential k-mer enrichment between conditions relative to the whole-genome background.

---

## Pipeline overview

```
01_find_threshold.R              # define region BEDs for each condition
        ↓
02_extract_sequences.sh          # extract FASTA sequences from GRCh38
        ↓
03_run_all_jellyfish.sh          # count k-mers (k=31) in hotspot regions + genome
  └── jellyfish_helper.sh        #   (called per-k by the above)
        ↓
04_kmer_differential_lowMem.py   # differential enrichment + statistical testing
        ↓
04_kmer_analysis.R               # visualization
```

---

## Scripts

### `01_find_threshold.R`

Sweeps an SV count threshold N to identify windows where `sv_count > N`, then measures base-pair overlap with the input region set. Computes precision, recall, and F1 at each threshold and selects the optimal N by F1. With `use_threshold = FALSE`, the threshold step is skipped and all input regions are written directly to the overlap BED.

**Inputs**
- `sv_windows_csv` — per-window SV counts (`*_varCount.csv`; columns: `chr`, `start`, `end`, `count`)
- `hotspot_path` — region set to analyze; accepts `.xlsx` or a 3-column BED/TSV

**Outputs** (written to `output_dir/`)
- `threshold_sweep.tsv` — precision / recall / F1 at every tested N
- `threshold_optimization.pdf` — sweep curves with optimal N annotated
- `high_sv_windows.bed` — windows passing the threshold
- `hotspot_high_sv_overlap.bed` — **primary output** passed to step 2

**Usage**

Edit the `main()` call at the bottom of the script. Run once per condition:

```r
main(sv_windows_csv = "...", hotspot_path = "sv_harbor_AllVar.bed",
     output_dir = "SVHarb_vs_AllGen_NoFilter", use_threshold = FALSE)

main(sv_windows_csv = "...", hotspot_path = "noTP_1_AllVar_inSVHarbor.bed",
     output_dir = "noTP_vs_AllGen_NoFilter", use_threshold = FALSE)
```

---

### `02_extract_sequences.sh`

Extracts FASTA sequences for the overlap regions produced in step 1. Indexes the genome with `samtools faidx` if needed, sorts the BED, and calls `bedtools getfasta`.

```bash
bash 02_extract_sequences.sh <genome.fa> <overlap.bed> <output.fa>
```

Run once per condition:

```bash
bash 02_extract_sequences.sh \
    GRCh38.fa \
    SVHarb_vs_AllGen_NoFilter/hotspot_high_sv_overlap.bed \
    SVHarb_vs_AllGen_NoFilter/hotspot_regions.fa

bash 02_extract_sequences.sh \
    GRCh38.fa \
    noTP_vs_AllGen_NoFilter/hotspot_high_sv_overlap.bed \
    noTP_vs_AllGen_NoFilter/hotspot_regions.fa
```

---

### `03_run_all_jellyfish.sh` + `jellyfish_helper.sh`

Counts canonical k-mers at k = 31 in both the hotspot FASTA and the full genome FASTA. `jellyfish_helper.sh` is called internally; it runs `jellyfish count` followed by `jellyfish dump` and writes a two-column TSV (`kmer<TAB>count`). Runs for a given k are skipped if the output TSV already exists.

```bash
bash 03_run_all_jellyfish.sh <hotspot_regions.fa> <genome.fa> <output_dir> [threads]
```

Run once per condition:

```bash
bash 03_run_all_jellyfish.sh \
    SVHarb_vs_AllGen_NoFilter/hotspot_regions.fa \
    GRCh38.fa \
    SVHarb_vs_AllGen_NoFilter \
    8

bash 03_run_all_jellyfish.sh \
    noTP_vs_AllGen_NoFilter/hotspot_regions.fa \
    GRCh38.fa \
    noTP_vs_AllGen_NoFilter \
    8
```

**Outputs** in each `output_dir/`:
- `hotspot_k31.tsv` — k-mer counts for the hotspot regions
- `genome_k31.tsv` — k-mer counts for the full genome

---

### `04_kmer_differential_lowMem.py`

For k = 31, tests whether any k-mer is differentially enriched in SVHarb hotspot regions vs. noTP hotspot regions, with each condition expressed relative to its own genome background.

**Statistical approach**

For each k-mer shared across all four count files, a 2×2 contingency table is constructed:

|                | this k-mer | all other k-mers |
|----------------|------------|-----------------|
| SVHarb hotspot | a          | b               |
| noTP hotspot   | c          | d               |

Chi-squared (no Yates correction) is used when all expected cell counts ≥ 5; otherwise falls back to Fisher's exact test (two-sided). Bonferroni correction is applied within each k independently.

Enrichment metadata (`log2_ratio`, `enrichment_or`) relative to the genome background is computed for each condition separately and included as output columns.

**Memory design:** count files are loaded as plain Python dicts rather than DataFrames, reducing peak memory by ~4–6× vs. a pandas-merge approach. Dicts are explicitly freed before the result DataFrame is built.

```bash
python3 04_kmer_differential_lowMem.py \
    --svharb_hotspot_prefix SVHarb_vs_AllGen_NoFilter/hotspot \
    --svharb_genome_prefix  SVHarb_vs_AllGen_NoFilter/genome  \
    --notp_hotspot_prefix   noTP_vs_AllGen_NoFilter/hotspot   \
    --notp_genome_prefix    noTP_vs_AllGen_NoFilter/genome    \
    --output_dir            Differential_SVHarb_vs_noTP       \
    --k_values 31
```

**Output** — `Differential_SVHarb_vs_noTP/kmer_differential_k31.tsv`, sorted by `pvalue_bonferroni` ascending. Columns:

| Column | Description |
|--------|-------------|
| `kmer` | k-mer sequence |
| `count_svharb`, `count_notp` | raw hotspot counts |
| `count_genome_sv`, `count_genome_notp` | raw genome counts |
| `log2_ratio_svharb`, `log2_ratio_notp` | log₂(hotspot freq / genome freq) per condition |
| `enrichment_or_svharb`, `enrichment_or_notp` | odds ratio vs. genome per condition |
| `delta_log2` | `log2_ratio_svharb − log2_ratio_notp` |
| `OR_differential` | odds ratio from the 2×2 differential table |
| `pvalue_raw`, `pvalue_bonferroni` | p-values |
| `test_used` | `chi2`, `fisher`, or `skipped` |

---

### `04_kmer_analysis.R`

Reads `kmer_differential_k31.tsv` and produces publication-ready plots. Does no computation — run `04_kmer_differential_lowMem.py` first.

Generates:
1. **Paired bar plot** — top N significant k-mers (Bonferroni < 0.05) showing log₂(hotspot/genome) for SVHarb and noTP, with significance annotations (`*` / `**` / `***`)
2. **Delta strip** — flanking panel showing Δlog₂ (noTP − SVHarb) on a diverging color scale

Configure `input_dir`, `output_dir`, `k_plot`, and `top_n_bars` in the user settings block at the top of the script before running.

---

## Dependencies

**R** (≥ 4.2): `tidyverse`, `readxl`, `valr`, `ggrepel`, `svglite`, `patchwork`

**Python** (≥ 3.10): `numpy`, `pandas`, `scipy`

**Shell**: `samtools`, `bedtools`, `jellyfish`

A Docker image with all shell dependencies pre-installed is available:

```
ztang301/all_dinumt:v1.1J
```

---

## Reference genome

All runs used GRCh38 no-alt analysis set:
`GCA_000001405.15_GRCh38_no_alt_analysis_set.fa`