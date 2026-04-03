#!/usr/bin/env bash
set -euo pipefail

# Run RTG vcfeval for SNV call sets across k-mer mappability regions.
#
# Usage:
#   bash 05_run_snv_compare_kmer.sh
#
# Optional environment overrides:
#   VERSION=v1.6
#   QUERY_BASE_DIR=/path/to/call_set_tumor_normal_revision
#   OUT_DIR=/path/to/output_dir
#   RELIABLE_DIR=/path/to/kmer_dir

VERSION="${VERSION:-v1.1.0}"

QUERY_BASE_DIR="${QUERY_BASE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/call_set_tumor_normal_revision}"
OUT_DIR="${OUT_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/snv/kmer}"
RELIABLE_DIR="${RELIABLE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/kmer}"

TRUTH_HG38="${TRUTH_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark_snvs.vcf.gz}"
TRUTH_CHM13="${TRUTH_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark_snvs.vcf.gz}"

SDF_HG38="${SDF_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/Benchmark_region/compare_fp/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf}"
SDF_CHM13="${SDF_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/GCA_009914755.4.chrNames.sdf}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

declare -A snv_files=(
    [clairS_GRCh38]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_clairS_GRCh38_snps.vcf.gz"
    [deepsomatic_CHM13]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_deepsomatic_CHM13_snps.vcf.gz"
    [deepsomatic_GRCh38]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_deepsomatic_GRCh38_snps.vcf.gz"
    [mutect2_CHM13]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_mutect2_CHM13_snps.vcf.gz"
    [mutect2_GRCh38]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_mutect2_GRCh38_snps.vcf.gz"
    [neusomatic_CHM13]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_neusomatic_CHM13_snps.vcf.gz"
    [neusomatic_GRCh38]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_neusomatic_GRCh38_snps.vcf.gz"
    [strelka2_CHM13]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_strelka2_CHM13_snps.vcf.gz"
    [strelka2_GRCh38]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_strelka2_GRCh38_snps.vcf.gz"
    [varscan2_CHM13]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_varscan2_CHM13_snps.vcf.gz"
    [varscan2_GRCh38]="${QUERY_BASE_DIR}/SNV/HapMap_Mixture_varscan2_GRCh38_snps.vcf.gz"
)

for variant in "${!snv_files[@]}"; do
    query="${snv_files[$variant]}"

    if [[ "${variant}" == *_CHM13 ]]; then
        sdf="${SDF_CHM13}"
        truth="${TRUTH_CHM13}"
        ref_name="chm13"
    else
        sdf="${SDF_HG38}"
        truth="${TRUTH_HG38}"
        ref_name="hg38"
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
    echo "Sample pair: HapMap_Mixture,${sample_name}"
    echo "========================================"

    for id in k24.umap k36.umap k50.umap k100.umap all; do
        rtg vcfeval -T 10 \
            --template "${sdf}" \
            --baseline "${truth}" \
            --calls "${query}" \
            --sample "HapMap_Mixture,${sample_name}" \
            --evaluation-regions "${RELIABLE_DIR}/${id}_unique_${ref_name}.bed" \
            --squash-ploidy \
            --all-records \
            --output "SNV_${variant}_kmer_${id}"
    done

    if [[ "${variant}" == *_CHM13 ]]; then
        rtg vcfeval -T 10 \
            --template "${sdf}" \
            --baseline "${truth}" \
            --calls "${query}" \
            --sample "HapMap_Mixture,${sample_name}" \
            --evaluation-regions "${RELIABLE_DIR}/all_unique_w_chm13only_chm13.bed" \
            --squash-ploidy \
            --all-records \
            --output "SNV_${variant}_kmer_all_w_chm13only"

        rtg vcfeval -T 10 \
            --template "${sdf}" \
            --baseline "${truth}" \
            --calls "${query}" \
            --sample "HapMap_Mixture,${sample_name}" \
            --evaluation-regions "${RELIABLE_DIR}/all_unique_wo_chm13only_chm13.bed" \
            --squash-ploidy \
            --all-records \
            --output "SNV_${variant}_kmer_all_wo_chm13only"
    fi
done

echo "All SNV k-mer evaluations completed."