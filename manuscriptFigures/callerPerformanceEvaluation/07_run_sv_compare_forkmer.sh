#!/usr/bin/env bash
#set -euo pipefail

# Run Truvari bench for SV call sets across k-mer mappability regions.
#
# Usage:
#   bash 07_run_sv_compare_kmer.sh
#
# Optional environment overrides:
#   VERSION=v1.6
#   QUERY_BASE_DIR=/path/to/call_set_tumor_normal_revision
#   OUT_DIR=/path/to/output_dir
#   RELIABLE_DIR=/path/to/kmer_dir

VERSION="${VERSION:-v1.1.0}"

QUERY_BASE_DIR="${QUERY_BASE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/call_set_tumor_normal_revision}"
OUT_DIR="${OUT_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/sv/kmer}"
RELIABLE_DIR="${RELIABLE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/kmer}"


TRUTH_HG38="${TRUTH_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/GRCh38/SMHTHAPMAP6_GRCh38_${VERSION}_somatic_benchmark_svs.vcf.gz}"
TRUTH_CHM13="${TRUTH_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/truth_set/v.1.1.0/CHM13/SMHTHAPMAP6_CHM13_${VERSION}_somatic_benchmark_svs.vcf.gz}"

REF_HG38="${REF_HG38:-/storage2/fs1/epigenome/Active/shared_smaht/References/SMAHT_References/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
REF_CHM13="${REF_CHM13:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/final_submission/data/reference/chm13/GCA_009914755.4.chrNames.fa}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

declare -A sv_files=(
    #[pbsv_CHM13]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_pbsv_CHM13_svs.vcf.gz"
    #[pbsv_GRCh38]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_pbsv_GRCh38_svs.vcf.gz"
    #[savana_CHM13]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_savana_CHM13_svs.vcf.gz"
    #[savana_GRCh38]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_savana_GRCh38_svs.vcf.gz"
    #[severus_CHM13]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_severus_CHM13_svs.vcf.gz"
    #[severus_GRCh38]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_severus_GRCh38_svs.vcf.gz"
    #[sniffles_CHM13]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_sniffles_CHM13_svs.vcf.gz"
    #[sniffles_GRCh38]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_sniffles_GRCh38_svs.vcf.gz"
    [svision_CHM13]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_svision_CHM13_svs.vcf.gz"
    #[svision_GRCh38]="${QUERY_BASE_DIR}/SV/HapMap_Mixture_svision_GRCh38_svs.vcf.gz"
)

for variant in "${!sv_files[@]}"; do
    query="${sv_files[$variant]}"

    if [[ "${variant}" == *_CHM13 ]]; then
        ref="${REF_CHM13}"
        truth="${TRUTH_CHM13}"
        ref_name="chm13"
    else
        ref="${REF_HG38}"
        truth="${TRUTH_HG38}"
        ref_name="hg38"
    fi

    if [[ "${variant}" == pbsv_* ]]; then
        sample_name="HapMap_Mixture_merged"

    elif [[ "${variant}" == sniffles_* ]]; then
        sample_name="SAMPLE"

    elif [[ "${variant}" == severus_CHM13 ]]; then
        sample_name="SMHTHAPMAP6-X-X-NN-B001-bcm-broad-washu-pbmm2_1.13.0_GRCh38.aligned.sorted.merged.reheader.realigned_chm13_dac"

    elif [[ "${variant}" == severus_GRCh38 ]]; then
        sample_name="SMHTHAPMAP6-X-X-NN-B001-bcm-broad-washu-pbmm2_1.13.0_GRCh38.aligned.sorted.merged"

    elif [[ "${variant}" == savana_CHM13 ]]; then
        sample_name="SMHTHAPMAP6-X-X-NN-B001-bcm-broad-washu-pbmm2_1.13.0_GRCh38.aligned.sorted.merged.reheader.realigned_chm13_dac"

    elif [[ "${variant}" == savana_GRCh38 ]]; then
        sample_name="SMHTHAPMAP6-X-X-NN-B001-bcm-broad-washu-pbmm2_1.13.0_GRCh38.aligned.sorted.merged"

    elif [[ "${variant}" == svision_CHM13 ]]; then
        # tumor sample = first sample in the VCF
        sample_name="SMHTHAPMAP6-X-X-NN-B001-bcm-broad-washu-pbmm2_1.13.0_GRCh38.aligned.sorted.merged.reheader.realigned_chm13_dac.bam"
    elif [[ "${variant}" == svision_GRCh38 ]]; then
        # tumor sample = first sample in the VCF
        sample_name="SMHTHAPMAP6-X-X-NN-B001-bcm-broad-washu-pbmm2_1.13.0_GRCh38.aligned.sorted.merged.bam"

    else
        echo "Unsupported variant key: ${variant}" >&2
        exit 1
    fi

    echo "========================================"
    echo "Processing ${variant}"
    echo "Query: ${query}"
    echo "Truth: ${truth}"
    echo "Sample pair: HapMap_Mixture,${sample_name}"
    echo "========================================"

    # Total performance across all benchmarkable regions
    truvari bench \
        -f "${ref}" \
        -b "${truth}" \
        -c "${query}" \
        --bSample "HapMap_Mixture" \
        --cSample "${sample_name}" \
        -s 50 \
        -S 50 \
        --typeignore \
        --pctseq 0 \
        --pick multi \
        --pctsize 0.7 \
        --includebed "${RELIABLE_DIR}/all_unique_${ref_name}.bed" \
        -o "SV_${variant}_total"

    # k-mer stratified performance
    for id in k24.umap k36.umap k50.umap k100.umap all; do
        truvari bench \
            -f "${ref}" \
            -b "${truth}" \
            -c "${query}" \
            --bSample "HapMap_Mixture" \
            --cSample "${sample_name}" \
            -s 50 \
            -S 50 \
            --typeignore \
            --pctseq 0 \
            --pick multi \
            --pctsize 0.7 \
            --includebed "${RELIABLE_DIR}/${id}_unique_${ref_name}.bed" \
            -o "SV_${variant}_kmer_${id}"
    done

    if [[ "${variant}" == *_CHM13 ]]; then
        truvari bench \
            -f "${ref}" \
            -b "${truth}" \
            -c "${query}" \
            --bSample "HapMap_Mixture" \
            --cSample "${sample_name}" \
            -s 50 \
            -S 50 \
            --typeignore \
            --pctseq 0 \
            --pick multi \
            --pctsize 0.7 \
            --includebed "${RELIABLE_DIR}/all_unique_w_chm13only_chm13.bed" \
            -o "SV_${variant}_kmer_all_w_chm13only"

        truvari bench \
            -f "${ref}" \
            -b "${truth}" \
            -c "${query}" \
            --bSample "HapMap_Mixture" \
            --cSample "${sample_name}" \
            -s 50 \
            -S 50 \
            --typeignore \
            --pctseq 0 \
            --pick multi \
            --pctsize 0.7 \
            --includebed "${RELIABLE_DIR}/all_unique_wo_chm13only_chm13.bed" \
            -o "SV_${variant}_kmer_all_wo_chm13only"
    fi
done

echo "All SV k-mer evaluations completed."