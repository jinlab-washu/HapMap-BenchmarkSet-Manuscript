#!/usr/bin/env bash
set -euo pipefail
#################################################################
# Split query-set body text by VAF and rebuild compressed/indexed VCFs.
#
# Usage:
#   bash 04_make_vaf_binned_queryset_vcfs.sh
#
# Optional environment overrides:
#   BASE_DIR=/path/to/call_set_tumor_normal_revision
#   REHEADER_DIR=/path/to/call_set_tumor_normal_revision/reheader
#   SNP_CODE_DIR=/path/to/code_from_zilan/rtg/snp
#   SV_CODE_DIR=/path/to/code_from_zilan/rtg/SV
#   TASK_IDS="1 2 3 4 5 6 7 8 9"
# Used Docker:cjyoon/ascp_triomix:v0.0.2a
# 
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu

BASE_DIR="${BASE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/call_set_tumor_normal_revision}"
REHEADER_DIR="${REHEADER_DIR:-${BASE_DIR}/reheader}"
CODE_DIR="${CODE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/code}"
TASK_IDS="${TASK_IDS:-1 2 3 4 5 6 7 8 9}"

declare -A snv_files=(
    [clairS_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_clairS_GRCh38_snps.vcf.gz"
    [deepsomatic_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_deepsomatic_CHM13_snps.vcf.gz"
    [deepsomatic_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_deepsomatic_GRCh38_snps.vcf.gz"
    [mutect2_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_mutect2_CHM13_snps.vcf.gz"
    [mutect2_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_mutect2_GRCh38_snps.vcf.gz"
    [neusomatic_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_neusomatic_CHM13_snps.vcf.gz"
    [neusomatic_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_neusomatic_GRCh38_snps.vcf.gz"
    [strelka2_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_strelka2_CHM13_snps.vcf.gz"
    [strelka2_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_strelka2_GRCh38_snps.vcf.gz"
    [varscan2_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_varscan2_CHM13_snps.vcf.gz"
    [varscan2_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_varscan2_GRCh38_snps.vcf.gz"
)

declare -A indel_files=(
    [deepsomatic_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_deepsomatic_CHM13_indels.vcf.gz"
    [deepsomatic_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_deepsomatic_GRCh38_indels.vcf.gz"
    [mutect2_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_mutect2_GRCh38_indels.vcf.gz"
    [neusomatic_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_neusomatic_CHM13_indels.vcf.gz"
    [neusomatic_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_neusomatic_GRCh38_indels.vcf.gz"
    [strelka2_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_strelka2_CHM13_indels.vcf.gz"
    [strelka2_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_strelka2_GRCh38_indels.vcf.gz"
    [varscan2_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_varscan2_CHM13_indels.vcf.gz"
    [varscan2_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_varscan2_GRCh38_indels.vcf.gz"
)

declare -A sv_files=(
    [pbsv_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_pbsv_CHM13_svs.vcf.gz"
    [pbsv_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_pbsv_GRCh38_svs.vcf.gz"
    [savana_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_savana_CHM13_svs.vcf.gz"
    [savana_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_savana_GRCh38_svs.vcf.gz"
    [severus_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_severus_CHM13_svs.vcf.gz"
    [severus_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_severus_GRCh38_svs.vcf.gz"
    [sniffles_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_sniffles_CHM13_svs.vcf.gz"
    [sniffles_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_sniffles_GRCh38_svs.vcf.gz"
    [svision_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_svision_CHM13_svs.vcf.gz"
)

run_small_variant_group() {
    local suffix="$1"
    local assoc_name="$2"
    declare -n files_ref="${assoc_name}"

    for key in "${!files_ref[@]}"; do
        echo "Processing ${suffix} for key: ${key}"

        Rscript "${CODE_DIR}/filter_AF_query_snvindel.R"             -f "${REHEADER_DIR}/${key}_${suffix}.txt"             -t 10             -n "${key}"             -c "${suffix}"             -o "${REHEADER_DIR}"

        for id in ${TASK_IDS}; do
            cat "${REHEADER_DIR}/${key}_${suffix}_header.txt"                 "${REHEADER_DIR}/${key}_${suffix}_${id}.txt"                 > "${REHEADER_DIR}/${key}_${suffix}_${id}.vcf"

            bgzip -c "${REHEADER_DIR}/${key}_${suffix}_${id}.vcf"                 > "${REHEADER_DIR}/${key}_${suffix}_${id}.vcf.gz"

            tabix -p vcf "${REHEADER_DIR}/${key}_${suffix}_${id}.vcf.gz"
        done
    done
}

run_sv_group() {
    for key in "${!sv_files[@]}"; do
        echo "Processing svs for key: ${key}"

        Rscript "${CODE_DIR}/filter_AF_sniffles.R"             -f "${REHEADER_DIR}/${key}_svs.txt"             -t 10             -n "${key}"             -c "svs"             -o "${REHEADER_DIR}"

        for id in ${TASK_IDS}; do
            cat "${REHEADER_DIR}/${key}_svs_header.txt"                 "${REHEADER_DIR}/${key}_svs_${id}.txt"                 > "${REHEADER_DIR}/${key}_svs_${id}.vcf"

            bgzip -c "${REHEADER_DIR}/${key}_svs_${id}.vcf"                 > "${REHEADER_DIR}/${key}_svs_${id}.vcf.gz"

            tabix -p vcf "${REHEADER_DIR}/${key}_svs_${id}.vcf.gz"
        done
    done
}

run_small_variant_group "snvs" "snv_files"
run_small_variant_group "indels" "indel_files"
run_sv_group
