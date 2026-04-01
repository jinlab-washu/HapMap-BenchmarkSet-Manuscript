#!/bin/bash
# 02_extract_sequences.sh
# Extract FASTA sequences for the hotspot/high-SV overlap regions
#
# Second step after you run /Github_Repos/code_backup_RIS/temp_code_local/202406_SMaHTmidYearMeeting/01_find_threshold.R
#
# Usage: bash 02_extract_sequences.sh <genome.fa> <overlap.bed> <output.fa>
# Requires: bedtools, samtools (for faidx index)

set -euo pipefail

GENOME_FA=${1:?  "arg1: genome FASTA (e.g. hg38.fa)"}
OVERLAP_BED=${2:?"arg2: overlap BED (output of 01_find_threshold.R)"}
OUTPUT_FA=${3:?  "arg3: output FASTA path"}

# Index genome if not already done
if [[ ! -f "${GENOME_FA}.fai" ]]; then
  echo "[index] Creating FASTA index for ${GENOME_FA} …"
  samtools faidx "$GENOME_FA"
fi

# Sort BED (bedtools getfasta requires sorted input matching FASTA order)
SORTED_BED=$(mktemp --suffix=.bed)
bedtools sort -i "$OVERLAP_BED" > "$SORTED_BED"

echo "[getfasta] Extracting sequences from ${OVERLAP_BED} …"
bedtools getfasta \
  -fi  "$GENOME_FA" \
  -bed "$SORTED_BED" \
  -fo  "$OUTPUT_FA"

rm -f "$SORTED_BED"

N_SEQS=$(grep -c "^>" "$OUTPUT_FA")
TOTAL_BP=$(grep -v "^>" "$OUTPUT_FA" | tr -d '\n' | wc -c)
echo "[done] ${N_SEQS} sequences  |  ${TOTAL_BP} bp  →  ${OUTPUT_FA}"