#!/usr/bin/env bash
set -euo pipefail

# Run Truvari bench for MEI call sets against truth sets.
#
# Usage:
#   bash 08_run_mei_compare.sh
#
# Optional environment overrides:
#   VERSION=v1.1.0
#   OUT_DIR=/path/to/output_dir
#   REF_HG38=/path/to/hg38.fa
#   REF_CHM13=/path/to/chm13.fa
#   RELIABLE_HG38=/path/to/hg38_reliable.bed
#   RELIABLE_CHM13=/path/to/chm13_reliable.bed
#
# Used Docker: meredith705/truvari:latest (or equivalent Truvari environment)
#
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

VERSION="${VERSION:-v1.1.0}"
OUT_DIR="${OUT_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/mei}"

# Truth sets
TRUTH_ALU_HG38="${TRUTH_ALU_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark_mei_alu.vcf.gz}"
TRUTH_L1_HG38="${TRUTH_L1_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark_mei_l1.vcf.gz}"
TRUTH_SVA_HG38="${TRUTH_SVA_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark_mei_sva.vcf.gz}"

TRUTH_ALU_CHM13="${TRUTH_ALU_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark_mei_alu.vcf.gz}"
TRUTH_L1_CHM13="${TRUTH_L1_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark_mei_l1.vcf.gz}"
TRUTH_SVA_CHM13="${TRUTH_SVA_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark_mei_sva.vcf.gz}"

# Query sets
QUERY_PALMER_ALU_HG38="${QUERY_PALMER_ALU_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_palmer/reformatted_sr0+/PALMER_ALU_hg38_reform_sorted.vcf.gz}"
QUERY_PALMER_L1_HG38="${QUERY_PALMER_L1_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_palmer/reformatted_sr0+/PALMER_LINE_hg38_reform_sorted.vcf.gz}"
QUERY_PALMER_SVA_HG38="${QUERY_PALMER_SVA_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_palmer/reformatted_sr0+/PALMER_SVA_hg38_reform_sorted.vcf.gz}"

QUERY_XTEA_ALU_HG38="${QUERY_XTEA_ALU_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_xTea/case_control_20250908/ALU_hg38/fn.vcf.gz}"
QUERY_XTEA_L1_HG38="${QUERY_XTEA_L1_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_xTea/case_control_20250908/LINE1_hg38/fn.vcf.gz}"
QUERY_XTEA_SVA_HG38="${QUERY_XTEA_SVA_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_xTea/case_control_20250908/SVA_hg38/fn.vcf.gz}"

QUERY_PALMER_ALU_CHM13="${QUERY_PALMER_ALU_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_palmer/reformatted_sr0+/PALMER_ALU_chm13_reform_sorted.vcf.gz}"
QUERY_PALMER_L1_CHM13="${QUERY_PALMER_L1_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_palmer/reformatted_sr0+/PALMER_LINE_chm13_reform_sorted.vcf.gz}"
QUERY_PALMER_SVA_CHM13="${QUERY_PALMER_SVA_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_palmer/reformatted_sr0+/PALMER_SVA_chm13_reform_sorted.vcf.gz}"

QUERY_XTEA_ALU_CHM13="${QUERY_XTEA_ALU_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_xTea/case_control_20250908/ALU_chm13/fn.vcf.gz}"
QUERY_XTEA_L1_CHM13="${QUERY_XTEA_L1_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_xTea/case_control_20250908/LINE1_chm13/fn.vcf.gz}"
QUERY_XTEA_SVA_CHM13="${QUERY_XTEA_SVA_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/zt_temp/reformatMEIcalls_20250311/reform_xTea/case_control_20250908/SVA_chm13/fn.vcf.gz}"
# References and benchmark regions
REF_HG38="${REF_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
REF_CHM13="${REF_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/data/reference/chm13/GCA_009914755.4.chrNames.fa}"

RELIABLE_HG38="${RELIABLE_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark.bed}"
RELIABLE_CHM13="${RELIABLE_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark.bed}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

declare -A mei_files=(
    [palmer_alu_GRCh38]="${QUERY_PALMER_ALU_HG38}"
    [palmer_l1_GRCh38]="${QUERY_PALMER_L1_HG38}"
    [palmer_sva_GRCh38]="${QUERY_PALMER_SVA_HG38}"
    [xtea_alu_GRCh38]="${QUERY_XTEA_ALU_HG38}"
    [xtea_l1_GRCh38]="${QUERY_XTEA_L1_HG38}"
    [xtea_sva_GRCh38]="${QUERY_XTEA_SVA_HG38}"

    [palmer_alu_CHM13]="${QUERY_PALMER_ALU_CHM13}"
    [palmer_l1_CHM13]="${QUERY_PALMER_L1_CHM13}"
    [palmer_sva_CHM13]="${QUERY_PALMER_SVA_CHM13}"
    [xtea_alu_CHM13]="${QUERY_XTEA_ALU_CHM13}"
    [xtea_l1_CHM13]="${QUERY_XTEA_L1_CHM13}"
    [xtea_sva_CHM13]="${QUERY_XTEA_SVA_CHM13}"
)

get_truth_vcf() {
    local variant="$1"

    case "${variant}" in
        *_alu_GRCh38)  echo "${TRUTH_ALU_HG38}" ;;
        *_l1_GRCh38)   echo "${TRUTH_L1_HG38}" ;;
        *_sva_GRCh38)  echo "${TRUTH_SVA_HG38}" ;;
        *_alu_CHM13)   echo "${TRUTH_ALU_CHM13}" ;;
        *_l1_CHM13)    echo "${TRUTH_L1_CHM13}" ;;
        *_sva_CHM13)   echo "${TRUTH_SVA_CHM13}" ;;
        *)
            echo "Unsupported MEI truth key: ${variant}" >&2
            exit 1
            ;;
    esac
}

get_reference_fasta() {
    local variant="$1"

    if [[ "${variant}" == *_CHM13 ]]; then
        echo "${REF_CHM13}"
    else
        echo "${REF_HG38}"
    fi
}

get_reliable_bed() {
    local variant="$1"

    if [[ "${variant}" == *_CHM13 ]]; then
        echo "${RELIABLE_CHM13}"
    else
        echo "${RELIABLE_HG38}"
    fi
}

get_query_sample_name() {
    local variant="$1"

    case "${variant}" in
        palmer_alu_GRCh38)  echo "PALMER_ALU_hg38" ;;
        palmer_l1_GRCh38)   echo "PALMER_LINE_hg38" ;;
        palmer_sva_GRCh38)  echo "PALMER_SVA_hg38" ;;
        xtea_alu_GRCh38)    echo "xTea_ALU_hg38_tumor" ;;
        xtea_l1_GRCh38)     echo "xTea_LINE1_hg38_tumor" ;;
        xtea_sva_GRCh38)    echo "xTea_SVA_hg38_tumor" ;;

        palmer_alu_CHM13)   echo "PALMER_ALU_chm13" ;;
        palmer_l1_CHM13)    echo "PALMER_LINE_chm13" ;;
        palmer_sva_CHM13)   echo "PALMER_SVA_chm13" ;;
        xtea_alu_CHM13)     echo "xTea_ALU_chm13_tumor" ;;
        xtea_l1_CHM13)      echo "xTea_LINE1_chm13_tumor" ;;
        xtea_sva_CHM13)     echo "xTea_SVA_chm13_tumor" ;;

        *)
            echo "Unsupported MEI query key: ${variant}" >&2
            exit 1
            ;;
    esac
}

for variant in "${!mei_files[@]}"; do
    query="${mei_files[$variant]}"
    truth="$(get_truth_vcf "${variant}")"
    ref="$(get_reference_fasta "${variant}")"
    reliable="$(get_reliable_bed "${variant}")"
    sample_name="$(get_query_sample_name "${variant}")"

    echo "========================================"
    echo "Processing ${variant}"
    echo "Query: ${query}"
    echo "Truth: ${truth}"
    echo "Reference: ${ref}"
    echo "Reliable regions: ${reliable}"
    echo "Sample pair: HapMap_Mixture,${sample_name}"
    echo "========================================"

    truvari bench \
        -f "${ref}" \
        -b "${truth}" \
        -c "${query}" \
        --bSample "HapMap_Mixture" \
        --cSample "${sample_name}" \
        --includebed "${reliable}" \
        -s 50 \
        -S 50 \
        --pick multi \
        --pctsize 0.1 \
        --pctseq 0.1 \
        --refdist 50 \
        -o "MEI_${variant}_total"
done

echo "All MEI evaluations completed."