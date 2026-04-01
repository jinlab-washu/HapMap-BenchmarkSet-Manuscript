#!/bin/bash
# jellyfish_helper.sh
# Count and dump k-mers from a single FASTA for a single k value.
#
# Usage: bash jellyfish_helper.sh <input.fa> <output_prefix> <k> [threads] [hash_size]
#
# Outputs:
#   <prefix>_k<k>.jf   – Jellyfish binary hash
#   <prefix>_k<k>.tsv  – tab-separated: kmer  count

set -euo pipefail

FASTA=${1:?     "arg1: input FASTA"}
PREFIX=${2:?    "arg2: output prefix (e.g. output/hotspot)"}
K=${3:?         "arg3: kmer size"}
THREADS=${4:-4}
HASH_SIZE=${5:-"1G"}   # increase to 4G for whole-genome runs

JF_FILE="${PREFIX}_k${K}.jf"
TSV_FILE="${PREFIX}_k${K}.tsv"

echo "[jellyfish] k=${K}  input=$(basename ${FASTA})  threads=${THREADS}"

# Count k-mers  (-C = canonical, counts forward + reverse complement together)
jellyfish count \
  -m  "$K" \
  -s  "$HASH_SIZE" \
  -t  "$THREADS" \
  -C  \
  -o  "$JF_FILE" \
  "$FASTA"

# Dump as tab-separated text (kmer  count)
jellyfish dump \
  -c \
  -t \
  "$JF_FILE" \
  > "$TSV_FILE"

N_KMERS=$(wc -l < "$TSV_FILE")
echo "[done]  ${N_KMERS} distinct k-mers  →  ${TSV_FILE}"