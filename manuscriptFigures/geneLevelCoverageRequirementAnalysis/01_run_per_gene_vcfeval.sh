#!/usr/bin/env bash
set -euo pipefail
#################################################################
# Run per-gene RTG vcfeval benchmarking in chunked LSF job arrays.
#
# This script:
#   1. Counts valid BED entries from a gene BED file
#   2. Splits the workload into chunks based on cluster array limits
#   3. Submits LSF job arrays using Dockerized RTG Tools
#   4. Runs vcfeval per gene for multiple callsets
#
# Usage:
#   bash 01_run_per_gene_vcfeval.sh
#
# Optional environment overrides:
#   TRUTH=/path/to/VAF=1%truth.vcf.gz 
#   CALL_DIR=/path/to/callsets
#   SDF=/path/to/reference.sdf
#   GENE_BED=/path/to/gene_regions.bed
#   RUN_DIR=/path/to/run_directory
#   WORK_ROOT=per_gene_vcfeval_runs_final
#   CONCURRENCY=300
#
# Used docker: blcdsdockerregistry/rtg-tools:3.12
#
# Author: Nahyun Kong
# Contact: nahyun@wustl.edu
#################################################################

export LSF_DOCKER_VOLUMES='/storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active /storage2/fs1/epigenome/Active/:/storage2/fs1/epigenome/Active/ /home/k.nahyun:/home/k.nahyun'

############################
# Input files and directories
############################

TRUTH="${TRUTH:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/Validation/truthset/one.vcf.gz}"
CALL_DIR="${CALL_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/data/geneset_analysis}"

BCM="${BCM:-${CALL_DIR}/mutect2_bcm.vcf.gz}"
NYGC="${NYGC:-${CALL_DIR}/mutect2_nygc.vcf.gz}"
BROAD="${BROAD:-${CALL_DIR}/mutect2_broad.vcf.gz}"
WASHU="${WASHU:-${CALL_DIR}/mutect2_washu.vcf.gz}"

SDF="${SDF:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/Benchmark_region/compare_fp/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf}"
GENE_BED="${GENE_BED:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/Validation/model_fitting_validation/per_gene_set/data/combined_gene.bed}"

############################
# Output directories
############################

RUN_DIR="${RUN_DIR:-/storage2/fs1/epigenome/Active/shared_smaht/TEST_RUN_FOLDER/TrueMutSet/nahyun_tmp/caller_evaluation/snv/per_gene/VAF1}"
WORK_ROOT="${WORK_ROOT:-per_gene_vcfeval_runs_final}"
LOG_DIR="${WORK_ROOT}/logs"

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"
mkdir -p "${LOG_DIR}"

# Use absolute path for Docker safety
WORK_ROOT="$(readlink -f "${WORK_ROOT}")"
export TRUTH SDF WASHU BCM NYGC BROAD GENE_BED WORK_ROOT

############################
# Determine job sizing
############################

N_GENES=$(awk 'NF>=4 && $1 !~ /^#/' "${GENE_BED}" | wc -l | tr -d ' ')
echo "Total genes: ${N_GENES}"

MAX_ARRAY=$(bparams -l 2>/dev/null | awk '/MAX_JOB_ARRAY_SIZE/ {print $3}' || true)
MAX_ARRAY="${MAX_ARRAY:-1000}"

if [[ -z "${MAX_ARRAY}" ]]; then
    MAX_ARRAY=1000
fi

CHUNK_SIZE=$(( MAX_ARRAY > 50 ? MAX_ARRAY - 50 : MAX_ARRAY ))
CONCURRENCY="${CONCURRENCY:-300}"

echo "Using chunk size: ${CHUNK_SIZE}"
echo "Cluster MAX_JOB_ARRAY_SIZE: ${MAX_ARRAY}"
echo "Concurrency cap: ${CONCURRENCY}"

############################
# Submit chunked LSF arrays
############################

start_idx=1
chunk_id=1

while (( start_idx <= N_GENES )); do
    end_idx=$(( start_idx + CHUNK_SIZE - 1 ))
    if (( end_idx > N_GENES )); then
        end_idx="${N_GENES}"
    fi

    chunk_size=$(( end_idx - start_idx + 1 ))

    echo "Submitting chunk ${chunk_id}: indices ${start_idx}-${end_idx}"

    bsub \
        -G compute-jin810-t3 \
        -q subscription \
        -sla jin810_t3 \
        -R 'rusage[mem=40GB]' \
        -M 40GB \
        -a 'docker(blcdsdockerregistry/rtg-tools:3.12)' \
        -J "vcfeval_per_gene_${chunk_id}[1-${chunk_size}]%${CONCURRENCY}" \
        -oo "${LOG_DIR}/chunk${chunk_id}.%I.out" \
        -eo "${LOG_DIR}/chunk${chunk_id}.%I.err" \
        /bin/bash -lc '
set -euo pipefail

offset=$(( '"${start_idx}"' - 1 ))
global_idx=$(( offset + LSB_JOBINDEX ))

line=$(awk -v tgt="$global_idx" '"'"'NF>=4 && $1 !~ /^#/ { if (++n==tgt) { print $1, $2, $3, $4; exit } }'"'"' "${GENE_BED}")

if [[ -z "$line" ]]; then
    echo "[warn] No BED line for index $global_idx" >&2
    exit 0
fi

set -- $line
chr=$1
startp=$2
endp=$3
gene=$4

tmp_bed=$(mktemp)
printf "%s\t%s\t%s\n" "$chr" "$startp" "$endp" > "$tmp_bed"

for center in washu BCM NYGC UW; do
    case "$center" in
        washu) calls="${WASHU}" ;;
        BCM)   calls="${BCM}" ;;
        NYGC)  calls="${NYGC}" ;;
        UW)    calls="${BROAD}" ;;
    esac

    outdir="${WORK_ROOT}/${gene}/${center}"
    rm -rf "$outdir"
    mkdir -p "$(dirname "$outdir")"

    rtg vcfeval \
        -f \
        -T 8 \
        --template "${SDF}" \
        --baseline "${TRUTH}" \
        --calls "$calls" \
        --sample HapMap_Mixture,TUMOR \
        --evaluation-regions "$tmp_bed" \
        --squash-ploidy \
        --all-records \
        --output "$outdir" \
        > "${outdir}.log" 2>&1 || true

    if [[ -s "${outdir}/summary.txt" ]]; then
        : > "${outdir}/.ok"
    else
        : > "${outdir}/.fail"
    fi
done

rm -f "$tmp_bed"
'
    start_idx=$(( end_idx + 1 ))
    chunk_id=$(( chunk_id + 1 ))
done

echo "Submitted all chunks."
