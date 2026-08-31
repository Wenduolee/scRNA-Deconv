#!/usr/bin/env bash
set -euo pipefail


# User Configurations & File Paths
REF="<path/to/reference_genome.fa>"
BAM="<path/to/aligned_sorted.bam>"
BARCODES="<path/to/barcodes.tsv>"
OUT_DIR="<path/to/output_directory>"

# Workflow Parameters
THREADS=24
CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"

mkdir -p "${OUT_DIR}"


# 1. Variant Calling (FreeBayes)

freebayes -f "${REF}" \
    --min-alternate-fraction 0.08 \
    --min-alternate-count 3 \
    --min-base-quality 20 \
    --min-mapping-quality 20 \
    --min-coverage 8 \
    "${BAM}" \
    > "${OUT_DIR}/raw_variants.vcf"


# 2. VCF Filtering & Indexing (bcftools / htslib)

bcftools view -i 'QUAL>20 & DP>10' "${OUT_DIR}/raw_variants.vcf" \
    | bgzip > "${OUT_DIR}/high_quality_snps.vcf.gz"

tabix -p vcf "${OUT_DIR}/high_quality_snps.vcf.gz"


# 3. Single-Cell Genotype Pileup (cellsnp-lite)

cellsnp-lite \
    -s "${BAM}" \
    -O "${OUT_DIR}/cellsnp_out" \
    -R "${OUT_DIR}/high_quality_snps.vcf.gz" \
    --barcodeFile "${BARCODES}" \
    --cellTAG CB \
    --UMItag UB \
    --minMAF 0.05 \
    --minCOUNT 5 \
    --minLEN 20 \
    --minMAPQ 10 \
    -p "${THREADS}" \
    --chrom "${CHROMS}" \
    --genotype \
    --gzip

# ==========================================
# 4. Donor Deconvolution (Vireo)
# ==========================================
vireo \
    -c "${OUT_DIR}/cellsnp_out" \
    -N "${N_DONORS}" \
    -o "${OUT_DIR}/vireo_out"
