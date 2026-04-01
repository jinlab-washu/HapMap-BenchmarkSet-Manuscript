#!/usr/bin/env bash
set -euo pipefail
#################################################################
# Summarize per-gene sensitivity from RTG vcfeval summary files.
#
# This script:
#   1. Reads gene names from a BED file
#   2. Looks for per-gene vcfeval summary.txt outputs
#   3. Extracts truth-set count and sensitivity for each center
#   4. Writes a combined TSV table
#
# Usage:
#   bash 02_collect_per_gene_sensitivity.sh
#
# Optional environment overrides:
#   WORK_ROOT=./per_gene_vcfeval_runs_final_totalsnv_include_nested
#   GENE_BED=/path/to/combined_gene.bed
#   OUT_TSV=/path/to/per_gene_sensitivity.tsv
#
# Expected input structure:
#   ${WORK_ROOT}/${gene}/${center}/summary.txt
#
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

WORK_ROOT="${WORK_ROOT:-././per_gene_vcfeval_runs_final}"
GENE_BED="${GENE_BED:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/Validation/model_fitting_validation/per_gene_set/data/combined_gene.bed}"
OUT_TSV="${OUT_TSV:-${WORK_ROOT}/per_gene_sensitivity.tsv}"

mkdir -p "${WORK_ROOT}"

parse_summary() {
    local summary_file="$1"
    awk '$1=="None" { tpb=$2; fn=$5; sens=$7; print tpb+fn, sens }' "${summary_file}"
}

echo -e "gene\tTruthset_number\tSensitivity_washu\tSensitivity_BCM\tSensitivity_NYGC\tSensitivity_BROAD" > "${OUT_TSV}"

awk 'NF>=4 && $1 !~ /^#/{print $4}' "${GENE_BED}" | while read -r gene; do
    truth="0"
    declare -A sens=()

    for center in washu BCM NYGC BROAD; do
        summary_file="${WORK_ROOT}/${gene}/${center}/summary.txt"

        if [[ -s "${summary_file}" ]]; then
            read -r truth sens["$center"] < <(parse_summary "${summary_file}")
        else
            sens["$center"]="NA"
        fi
    done

    echo -e "${gene}\t${truth}\t${sens[washu]}\t${sens[BCM]}\t${sens[NYGC]}\t${sens[BROAD]}" >> "${OUT_TSV}"
done

echo "[done] ${OUT_TSV}"
