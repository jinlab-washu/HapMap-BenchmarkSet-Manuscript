#!/bin/bash
# 03_run_all_jellyfish.sh
# Run Jellyfish for k=2..6 on both the hotspot regions and the full genome.
#
# Usage: bash 03_run_all_jellyfish.sh \
#            <hotspot_regions.fa> \
#            <genome.fa> \
#            <output_dir> \
#            [threads]
#
# Produces in <output_dir>/:
#   hotspot_k{2..6}.tsv   – hotspot region k-mer counts
#   genome_k{2..6}.tsv    – whole-genome k-mer counts

set -euo pipefail

HOTSPOT_FA=${1:?  "arg1: hotspot regions FASTA (from 02_extract_sequences.sh)"}
GENOME_FA=${2:?   "arg2: full genome FASTA (hg38.fa)"}
OUTPUT_DIR=${3:?  "arg3: output directory"}
THREADS=${4:-4}

mkdir -p "$OUTPUT_DIR"
SCRIPT_DIR="$(dirname "$0")"

K_VALUES=($(seq 2 31)) # prev: (2 3 4 5 6)

echo "===== Hotspot regions ====="
for K in "${K_VALUES[@]}"; do
  ## non-interactive mode
  TSV_FILE="${OUTPUT_DIR}/hotspot_k${K}.tsv"

  if [[ -f "$TSV_FILE" ]]; then
    echo "[skip] ${TSV_FILE} already exists — skipping hotspot k=${K}"
    continue
  fi

  echo "[process] hotspot k=${K}"
  bsub \
    -J "jellyfish_hotspot_k${K}" \
    -G compute-jin810 \
    -q general \
    -R 'rusage[mem=100GB]' \
    -n "$THREADS" \
    -a 'docker(ztang301/all_dinumt:v1.1J)' \
    bash "${SCRIPT_DIR}/jellyfish_helper.sh" \
      "$HOTSPOT_FA" \
      "${OUTPUT_DIR}/hotspot" \
      "$K" \
      "$THREADS" \
      "10G"
  
  ## Interactive mode
  # bash "$(dirname "$0")/jellyfish_helper.sh" \
  #   "$HOTSPOT_FA" \
  #   "${OUTPUT_DIR}/hotspot" \
  #   "$K" \
  #   "$THREADS" \
  #   "10G"   # 500M used for k=2-6, increase due to larger ks
done

echo ""
echo "===== Full genome ====="
for K in "${K_VALUES[@]}"; do
  ## non-interactive mode
  TSV_FILE="${OUTPUT_DIR}/genome_k${K}.tsv"

  if [[ -f "$TSV_FILE" ]]; then
    echo "[skip] ${TSV_FILE} already exists — skipping genome k=${K}"
    continue
  fi

  echo "[process] genome k=${K}"
  bsub \
    -J "jellyfish_genome_k${K}" \
    -G compute-jin810 \
    -q general \
    -R 'rusage[mem=150GB]' \
    -n "$THREADS" \
    -a 'docker(ztang301/all_dinumt:v1.1J)' \
    bash "${SCRIPT_DIR}/jellyfish_helper.sh" \
      "$GENOME_FA" \
      "${OUTPUT_DIR}/genome" \
      "$K" \
      "$THREADS" \
      "10G"
  
  ## interactive mode
  # bash "$(dirname "$0")/jellyfish_helper.sh" \
  #   "$GENOME_FA" \
  #   "${OUTPUT_DIR}/genome" \
  #   "$K" \
  #   "$THREADS" \
  #   "10G"
done

echo ""
echo "All Jellyfish runs complete. TSVs in: ${OUTPUT_DIR}"
ls "${OUTPUT_DIR}"/*.tsv