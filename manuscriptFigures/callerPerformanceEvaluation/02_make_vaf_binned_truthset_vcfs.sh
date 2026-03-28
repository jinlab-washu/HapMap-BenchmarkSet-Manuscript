#!/usr/bin/env bash
set -euo pipefail
#################################################################
# Split truth-set text files by VAF and rebuild compressed/indexed VCFs.
#
# Usage:
#   bash 02_make_vaf_binned_truthset_vcfs.sh
#
# Optional environment overrides:
#   VERSION=v1.1.0
#   CODE_DIR=/path/to/code_from_zilan/rtg/truth_set
#   INPUT_DIR=/path/to/data/truth_set/reheader
#   OUTPUT_DIR=/path/to/data/truth_set/reheader
#   WORK_DIR=/path/to/data/truth_set
#   REFERENCES="GRCh38 CHM13"
#   TASK_IDS="1 2 3 4 5 6 7 8 9"
#
# Used Docker:cjyoon/ascp_triomix:v0.0.2a
# 
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

VERSION="${VERSION:-v1.1.0}"
CODE_DIR="${CODE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/code/}"
INPUT_DIR="${INPUT_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/truth_set/reheader}"
OUTPUT_DIR="${OUTPUT_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/truth_set/reheader}"
WORK_DIR="${WORK_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/truth_set}"
REFERENCES="${REFERENCES:-GRCh38 CHM13}"
TASK_IDS="${TASK_IDS:-1 2 3 4 5 6 7 8 9}"

cd "${WORK_DIR}"

for reference in ${REFERENCES}; do
    echo "=============================="
    echo "Processing reference: ${reference}"
    echo "=============================="

    for variant_class in snvs indels svs; do
        echo "[${reference}] Running filter_AF_truth.R for ${variant_class}"

        Rscript "${CODE_DIR}/filter_AF_truth.R" \
            -f "SMHTHAPMAP6_${reference}_${VERSION}_somatic_benchmark_${variant_class}.txt" \
            -t 10 \
            -r "${reference}" \
            -c "${variant_class}" \
            -v "${VERSION}" \
            -o "${OUTPUT_DIR}"

        header_file="SMHTHAPMAP6_${reference}_${VERSION}_somatic_benchmark_${variant_class}_header.txt"

        for id in ${TASK_IDS}; do
            echo "[${reference}] ${variant_class}: task ID ${id}"

            txt_file="${INPUT_DIR}/truth_set_${VERSION}_${variant_class}_${reference}_${id}.txt"
            vcf_file="${OUTPUT_DIR}/truth_set_${VERSION}_${variant_class}_${reference}_${id}.vcf"
            vcfgz_file="${vcf_file}.gz"

            cat "${header_file}" "${txt_file}" > "${vcf_file}"
            bgzip -c "${vcf_file}" > "${vcfgz_file}"
            tabix -p vcf "${vcfgz_file}"
        done
    done
done
