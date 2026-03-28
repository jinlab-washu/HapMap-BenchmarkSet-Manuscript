#!/usr/bin/env bash
set -euo pipefail
#################################################################
# Split query-set VCFs into header and body text files.
#
# Usage:
#   bash 03_split_queryset_vcfs.sh
#
# Optional environment overrides:
#   BASE_DIR=/path/to/call_set_tumor_normal_revision
#   REHEADER_DIR=/path/to/call_set_tumor_normal_revision/reheader
# Used Docker:mgibio/bcftools-cwl:1.12
# 
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

BASE_DIR="${BASE_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/call_set_tumor_normal_revision}"
REHEADER_DIR="${REHEADER_DIR:-${BASE_DIR}/reheader}"

mkdir -p "${REHEADER_DIR}"

declare -A snv_files=(
 #   [clairS_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_clairS_GRCh38_snps.vcf.gz"
 #   [deepsomatic_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_deepsomatic_CHM13_snps.vcf.gz"
 #   [deepsomatic_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_deepsomatic_GRCh38_snps.vcf.gz"
 #   [mutect2_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_mutect2_CHM13_snps.vcf.gz"
 #   [mutect2_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_mutect2_GRCh38_snps.vcf.gz"
    [neusomatic_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_neusomatic_CHM13_snps.vcf.gz"
    [neusomatic_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_neusomatic_GRCh38_snps.vcf.gz"
 #   [strelka2_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_strelka2_CHM13_snps.vcf.gz"
 #   [strelka2_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_strelka2_GRCh38_snps.vcf.gz"
 #   [varscan2_CHM13]="${BASE_DIR}/SNV/HapMap_Mixture_varscan2_CHM13_snps.vcf.gz"
 #   [varscan2_GRCh38]="${BASE_DIR}/SNV/HapMap_Mixture_varscan2_GRCh38_snps.vcf.gz"
)

declare -A indel_files=(
#    [deepsomatic_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_deepsomatic_CHM13_indels.vcf.gz"
#    [deepsomatic_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_deepsomatic_GRCh38_indels.vcf.gz"
#    [mutect2_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_mutect2_GRCh38_indels.vcf.gz"
    [mutect2_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_mutect2_CHM13_indels.vcf.gz"
#    [neusomatic_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_neusomatic_CHM13_indels.vcf.gz"
#    [neusomatic_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_neusomatic_GRCh38_indels.vcf.gz"
#    [strelka2_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_strelka2_CHM13_indels.vcf.gz"
#    [strelka2_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_strelka2_GRCh38_indels.vcf.gz"
#    [varscan2_CHM13]="${BASE_DIR}/INDEL/HapMap_Mixture_varscan2_CHM13_indels.vcf.gz"
#    [varscan2_GRCh38]="${BASE_DIR}/INDEL/HapMap_Mixture_varscan2_GRCh38_indels.vcf.gz"
)

declare -A sv_files=(
#    [pbsv_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_pbsv_CHM13_svs.vcf.gz"
#    [pbsv_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_pbsv_GRCh38_svs.vcf.gz"
#    [savana_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_savana_CHM13_svs.vcf.gz"
#    [savana_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_savana_GRCh38_svs.vcf.gz"
#    [severus_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_severus_CHM13_svs.vcf.gz"
#    [severus_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_severus_GRCh38_svs.vcf.gz"
#    [sniffles_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_sniffles_CHM13_svs.vcf.gz"
#    [sniffles_GRCh38]="${BASE_DIR}/SV/HapMap_Mixture_sniffles_GRCh38_svs.vcf.gz"
#    [svision_CHM13]="${BASE_DIR}/SV/HapMap_Mixture_svision_CHM13_svs.vcf.gz"
)

process_group() {
    local suffix="$1"
    local assoc_name="$2"
    declare -n files_ref="${assoc_name}"

    for key in "${!files_ref[@]}"; do
        local file="${files_ref[$key]}"
        local header_file="${REHEADER_DIR}/${key}_${suffix}_header.txt"
        local body_file="${REHEADER_DIR}/${key}_${suffix}.txt"

        echo "Processing ${suffix}: ${key}"
        bcftools view --threads 10 -h "${file}" > "${header_file}"
        bcftools view --threads 10 -H "${file}" > "${body_file}"
    done
}

#process_group "snvs" "snv_files"
process_group "indels" "indel_files"
#process_group "svs" "sv_files"
