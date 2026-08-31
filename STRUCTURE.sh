#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Pipeline: Biogeographic Ancestry Analysis & STRUCTURE 
# STRUCTURE & Harvester Parameters
BURNIN=10000                 # 10,000 burn-in iterations
NUMREPS=10000                # 10,000 MCMC steps
K_MIN=2                      # Tested K range: 2 to 6
K_MAX=6
RUNS_PER_K=5                 # Number of replicates per K (standard for Evanno Delta K calculation)
STRUCTURE_BIN="structure"    # Path to structure binary
HARVESTER_BIN="structureHarvester.py"
# ==============================================================================


# Configurable Parameters
SAMPLE_VCF="<path/to/deconvolved_donor_genotypes.vcf.gz>"   # Deconvolved VCF from Vireo
KG_REF_VCF="<path/to/1000G_Phase3_All_Chr.vcf.gz>"          # 1000 Genomes Phase 3 VCF
OUT_DIR="<path/to/ancestry_output_dir>"
PLINK_BIN="plink"                                            # Path to plink executable

# Filtering Thresholds
LD_WINDOW_KB=1000
LD_STEP=1
LD_R2=0.2
HWE_PVAL=0.05
DELTA_MAF=0.3

mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# Step 1: Standardize Sample VCF & Retain Biallelic SNPs
# ------------------------------------------------------------------------------

# Set canonical variant ID (CHROM:POS:REF:ALT) to resolve missing/dot IDs
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' \
    -Oz -o "${OUT_DIR}/sample_annotated.vcf.gz" \
    "${SAMPLE_VCF}"

# Filter for biallelic SNPs only
bcftools view -m2 -M2 -v snps \
    -Oz -o "${OUT_DIR}/sample_biallelic.vcf.gz" \
    "${OUT_DIR}/sample_annotated.vcf.gz"

bcftools index -f "${OUT_DIR}/sample_biallelic.vcf.gz"

# Convert to PLINK binary format
${PLINK_BIN} --vcf "${OUT_DIR}/sample_biallelic.vcf.gz" \
    --keep-allele-order \
    --make-bed \
    --out "${OUT_DIR}/sample_clean"

# ------------------------------------------------------------------------------
# Step 2: Intersect Sample SNPs with 1000 Genomes Reference
# ------------------------------------------------------------------------------

# Extract genomic coordinates (CHROM:POS)
bcftools query -f '%CHROM\t%POS\n' "${OUT_DIR}/sample_biallelic.vcf.gz" \
    > "${OUT_DIR}/sample_snps_coord.txt"

# Subset 1000 Genomes VCF by target coordinates
bcftools view -Oz -R "${OUT_DIR}/sample_snps_coord.txt" "${KG_REF_VCF}" \
    -o "${OUT_DIR}/1kg_intersected.vcf.gz"

# Convert intersected 1000G data to PLINK format
${PLINK_BIN} --vcf "${OUT_DIR}/1kg_intersected.vcf.gz" \
    --vcf-half-call h \
    --keep-allele-order \
    --make-bed \
    --out "${OUT_DIR}/1kg_clean"

# ------------------------------------------------------------------------------
# Step 3: Linkage Disequilibrium (LD) & Hardy-Weinberg Equilibrium (HWE) Pruning
# ------------------------------------------------------------------------------
# LD pruning across 1000G reference loci
${PLINK_BIN} --bfile "${OUT_DIR}/1kg_clean" \
    --indep-pairwise "${LD_WINDOW_KB}kb" "${LD_STEP}" "${LD_R2}" \
    --out "${OUT_DIR}/1kg_ld"

# HWE test on LD-pruned loci
${PLINK_BIN} --bfile "${OUT_DIR}/1kg_clean" \
    --extract "${OUT_DIR}/1kg_ld.prune.in" \
    --hardy \
    --out "${OUT_DIR}/1kg_hwe"

# Filter loci with HWE p-value > 0.05
awk -v pval="${HWE_PVAL}" '$9 > pval {print $2}' "${OUT_DIR}/1kg_hwe.hwe" \
    | awk '!seen[$0]++' > "${OUT_DIR}/hwe_pass_snps.txt"

# ------------------------------------------------------------------------------
# Step 4: Extract Shared AISNPs & Merge Datasets
# ------------------------------------------------------------------------------
# Filter both datasets to high-confidence pruned AISNPs
${PLINK_BIN} --bfile "${OUT_DIR}/1kg_clean" \
    --extract "${OUT_DIR}/hwe_pass_snps.txt" \
    --make-bed \
    --out "${OUT_DIR}/1kg_aisnps"

${PLINK_BIN} --bfile "${OUT_DIR}/sample_clean" \
    --extract "${OUT_DIR}/hwe_pass_snps.txt" \
    --make-bed \
    --out "${OUT_DIR}/sample_aisnps"

# Merge sample dataset with 1000 Genomes reference
${PLINK_BIN} --bfile "${OUT_DIR}/sample_aisnps" \
    --bmerge "${OUT_DIR}/1kg_aisnps" \
    --keep-allele-order \
    --make-bed \
    --out "${OUT_DIR}/combined_aisnps"

# ------------------------------------------------------------------------------
# Step 5: Export STRUCTURE Input File
# ------------------------------------------------------------------------------
# Export to STRUCTURE format (.strct_in)
${PLINK_BIN} --bfile "${OUT_DIR}/combined_aisnps" \
    --recode structure \
    --out "${OUT_DIR}/structure_input"


# ------------------------------------------------------------------------------
# Step 6: Automated STRUCTURE Execution (K = 2 to 6, Admixture Model)
# ------------------------------------------------------------------------------

STR_DIR="${OUT_DIR}/structure_runs"
mkdir -p "${STR_DIR}"

# Dynamically calculate sample count (NUMINDS) and marker count (NUMLOCI) from PLINK bim/fam
NUMINDS=$(wc -l < "${OUT_DIR}/combined_aisnps.fam")
NUMLOCI=$(wc -l < "${OUT_DIR}/combined_aisnps.bim")

# Generate mainparams configuration file
cat << EOF > "${STR_DIR}/mainparams"
#define OUTFILE ${STR_DIR}/str_out
#define INFILE ${OUT_DIR}/structure_input.strct_in
#define NUMINDS ${NUMINDS}
#define NUMLOCI ${NUMLOCI}
#define PLOIDY 2
#define MISSING -9
#define ONEROWPERIND 0
#define LABEL 1
#define POPDATA 0
#define POPFLAG 0
#define LOCDATA 0
#define PHENOTYPE 0
#define MARKERNAMES 1
#define MAPDISTANCES 0
#define MAXPOPS 2
#define BURNIN ${BURNIN}
#define NUMREPS ${NUMREPS}
EOF

# Generate extraparams configuration file (Admixture Model)
cat << EOF > "${STR_DIR}/extraparams"
#define NOADMIX 0
#define LINKAGE 0
#define USEPOPINFO 0
#define LOCPRIOR 0
#define FREQSCORR 1
#define ONEFST 0
#define INFERALPHA 1
#define POPALPHAS 0
#define ALPHA 1.0
#define INFERLAMBDA 0
#define POPSPECIFICLAMBDA 0
#define LAMBDA 1.0
EOF

# Execute STRUCTURE iterations across K = 2..6 and multiple replicates
for ((k=K_MIN; k<=K_MAX; k++)); do
    for ((r=1; r<=RUNS_PER_K; r++)); do
        echo "Running STRUCTURE: K=${k}, Iteration=${r}/${RUNS_PER_K}..."
        ${STRUCTURE_BIN} \
            -m "${STR_DIR}/mainparams" \
            -e "${STR_DIR}/extraparams" \
            -K "${k}" \
            -o "${STR_DIR}/k${k}_run${r}_out"
    done
done

# ------------------------------------------------------------------------------
# Step 7: Optimal K Determination (Structure Harvester / Evanno Method)
# ------------------------------------------------------------------------------

HARVESTER_OUT="${OUT_DIR}/structure_harvester_out"
mkdir -p "${HARVESTER_OUT}"

# Run Structure Harvester CLI to parse log-likelihood and generate Delta K plots
${HARVESTER_BIN} \
    --dir="${STR_DIR}" \
    --out="${HARVESTER_OUT}" \
    --evanno \
    --clumpp

echo "=== All Pipeline Steps Finished Successfully ==="
echo "STRUCTURE results saved to: ${STR_DIR}"
echo "Optimal K evaluation results saved to: ${HARVESTER_OUT}"
