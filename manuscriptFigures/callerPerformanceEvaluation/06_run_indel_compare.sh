#!/usr/bin/env bash
set -euo pipefail

# Run RTG vcfeval for INDEL call sets against truth sets.
#
# Usage:
#   bash 06_run_indel_compare.sh
#
# Optional environment overrides:
#   VERSION=v1.1.0
#   QUERY_BASE_DIR=/path/to/call_set_tumor_normal_revision
#   QUERY_DIR=/path/to/call_set_tumor_normal_revision/reheader
#   TRUTH_DIR=/path/to/truth_set/reheader
#   OUT_DIR=/path/to/output_dir
#   RELIABLE_HG38=/path/to/hg38_reliable.bed
#   RELIABLE_CHM13=/path/to/chm13_reliable.bed
# Used Docker:blcdsdockerregistry/rtg-tools:3.12
# 
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

VERSION="${VERSION:-v1.1.0}"

TRUTH_DIR="${TRUTH_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/truth_set/reheader}"
QUERY_BASE_DIR="${QUERY_BASE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/call_set_tumor_normal_revision}"
QUERY_DIR="${QUERY_DIR:-${QUERY_BASE_DIR}/reheader}"
OUT_DIR="${OUT_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/indels}"

TRUTH_HG38="${TRUTH_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark_indels.vcf.gz}"
TRUTH_CHM13="${TRUTH_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark_indels.vcf.gz}"

REF_HG38="${REF_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
SDF_HG38="${SDF_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/Benchmark_region/compare_fp/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf}"
RELIABLE_HG38="${RELIABLE_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark.bed}"

REF_CHM13="${REF_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/data/reference/chm13/GCA_009914755.4.chrNames.fa}"
SDF_CHM13="${SDF_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/GCA_009914755.4.chrNames.sdf}"
RELIABLE_CHM13="${RELIABLE_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark.bed}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

declare -A indel_files=(
    [deepsomatic_CHM13]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_deepsomatic_CHM13_indels.vcf.gz"
    [deepsomatic_GRCh38]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_deepsomatic_GRCh38_indels.vcf.gz"
    [mutect2_CHM13]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_mutect2_CHM13_indels.vcf.gz"
    [mutect2_GRCh38]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_mutect2_GRCh38_indels.vcf.gz"
    [neusomatic_CHM13]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_neusomatic_CHM13_indels.vcf.gz"
    [neusomatic_GRCh38]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_neusomatic_GRCh38_indels.vcf.gz"
    [strelka2_CHM13]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_strelka2_CHM13_indels.vcf.gz"
    [strelka2_GRCh38]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_strelka2_GRCh38_indels.vcf.gz"
    [varscan2_CHM13]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_varscan2_CHM13_indels.vcf.gz"
    [varscan2_GRCh38]="${QUERY_BASE_DIR}/INDEL/HapMap_Mixture_varscan2_GRCh38_indels.vcf.gz"
)

for variant in "${!indel_files[@]}"; do
    query="${indel_files[$variant]}"

    if [[ "${variant}" == *_CHM13 ]]; then
        sdf="${SDF_CHM13}"
        truth="${TRUTH_CHM13}"
        reliable="${RELIABLE_CHM13}"
        truth_af_prefix="${TRUTH_DIR}/truth_set_${VERSION}_indels_CHM13"
    else
        sdf="${SDF_HG38}"
        truth="${TRUTH_HG38}"
        reliable="${RELIABLE_HG38}"
        truth_af_prefix="${TRUTH_DIR}/truth_set_${VERSION}_indels_GRCh38"
    fi

    if [[ "${variant}" == mutect2_* ]] || [[ "${variant}" == deepsomatic_* ]]; then
        sample_name="Hapmap_Mixture"

    elif [[ "${variant}" == strelka2_* ]] || [[ "${variant}" == varscan2_* ]]; then
        sample_name="TUMOR"

    elif [[ "${variant}" == neusomatic_* ]] || [[ "${variant}" == clairS_* ]]; then
        sample_name="SAMPLE"

    else
        echo "Unknown sample type for ${variant}"
        exit 1
    fi

    echo "========================================"
    echo "Processing ${variant}"
    echo "Query: ${query}"
    echo "Truth: ${truth}"
    echo "Reliable regions: ${reliable}"
    echo "Sample pair: HapMap_Mixture,${sample_name}"
    echo "========================================"

    # Total performance
    rtg vcfeval -T 10 \
        --template "${sdf}" \
        --baseline "${truth}" \
        --calls "${query}" \
        --sample "HapMap_Mixture,${sample_name}" \
        --evaluation-regions "${reliable}" \
        --squash-ploidy \
        --output "INDEL_${variant}_total"

    # Sensitivity: truth AF bins vs full query set
    for id in {1..9}; do
        rtg vcfeval -T 10 \
            --template "${sdf}" \
            --baseline "${truth_af_prefix}_${id}.vcf.gz" \
            --calls "${query}" \
            --sample "HapMap_Mixture,${sample_name}" \
            --evaluation-regions "${reliable}" \
            --squash-ploidy \
            --output "INDEL_${variant}_af_${id}_sensitivity"
    done

    # Precision: full truth set vs query AF bins
    for id in {1..9}; do
        rtg vcfeval -T 10 \
            --template "${sdf}" \
            --baseline "${truth}" \
            --calls "${QUERY_DIR}/${variant}_indels_${id}.vcf.gz" \
            --sample "HapMap_Mixture,${sample_name}" \
            --evaluation-regions "${reliable}" \
            --squash-ploidy \
            --output "INDEL_${variant}_af_${id}_precision"
    done
done

echo "All INDEL evaluations completed."
