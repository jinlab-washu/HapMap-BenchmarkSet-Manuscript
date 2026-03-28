#!/usr/bin/env bash
set -euo pipefail
#################################################################
# Split benchmark VCFs into header and body text files for each reference.
#
# Usage:
#   bash 01_split_truthset_vcfs.sh
#
# Optional environment overrides:
#   VERSION=v1.1.0
#   BASE_TRUTH_DIR=/path/to/final_submission/truth_set
#   WORK_DIR=/path/to/data/truth_set
#   REFERENCES="GRCh38 CHM13"
# Used Docker:mgibio/bcftools-cwl:1.12
# 
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

VERSION="${VERSION:-v1.1.0}"
BASE_TRUTH_DIR="${BASE_TRUTH_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set}"
WORK_DIR="${WORK_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/truth_set}"
REFERENCES="${REFERENCES:-GRCh38 CHM13}"

cd "${WORK_DIR}"

for reference in ${REFERENCES}; do
    echo "========================================"
    echo "Splitting benchmark VCFs for ${reference}"
    echo "========================================"

    truth_dir="${BASE_TRUTH_DIR}/v.1.1.0/${reference}"

    declare -A vcf_paths=(
        [snvs]="${truth_dir}/SMHTHAPMAP6_${reference}_${VERSION}_somatic_benchmark_snvs.vcf.gz"
        [indels]="${truth_dir}/SMHTHAPMAP6_${reference}_${VERSION}_somatic_benchmark_indels.vcf.gz"
        [svs]="${truth_dir}/SMHTHAPMAP6_${reference}_${VERSION}_somatic_benchmark_svs.vcf.gz"
    )

    for variant_class in snvs indels svs; do
        input_vcf="${vcf_paths[$variant_class]}"
        out_prefix="SMHTHAPMAP6_${reference}_${VERSION}_somatic_benchmark_${variant_class}"

        echo "[${reference}] ${variant_class}: ${input_vcf}"

        bcftools view --threads 10 -h "${input_vcf}" > "${out_prefix}_header.txt"
        bcftools view --threads 10 -H "${input_vcf}" > "${out_prefix}.txt"
    done
done
