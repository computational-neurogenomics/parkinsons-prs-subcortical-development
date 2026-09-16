#!/bin/bash
# ==============================================================================
# PIPELINE: SBayesRC POLYGENIC RISK SCORE CALCULATION (GCTB + PLINK2)
# ==============================================================================

# User-defined paths and parameters
GWAS_SUMSTATS="test.ma"
LDM_EIGEN="ldm"
ANNOT_FILE="annot.txt"
QC_VARIANTS="/QCed_variants_MAF01_CR90_RSQ60.list"
IMPUTED_VCF_DIR="/ImputedFiles"
OUTPUT_DIR="PRScalc"
THREADS=4
MEM_MB=35000

mkdir -p ${OUTPUT_DIR}

# ------------------------------------------------------------------------------
# Step 1: Summary Statistics Imputation using GCTB
# ------------------------------------------------------------------------------
echo "Step 1: Imputing summary statistics with GCTB..."
gctb --ldm-eigen ${LDM_EIGEN} \
     --gwas-summary ${GWAS_SUMSTATS} \
     --impute-summary \
     --out test \
     --thread ${THREADS}

# ------------------------------------------------------------------------------
# Step 2: SBayesRC Posterior Effect Size Estimation with Annotations
# ------------------------------------------------------------------------------
echo "Step 2: Estimating posterior effect sizes via SBayesRC..."
gctb --ldm-eigen ${LDM_EIGEN} \
     --gwas-summary test.imputed.ma \
     --sbayes RC \
     --annot ${ANNOT_FILE} \
     --out trait_fullLDM \
     --thread ${THREADS}

# ------------------------------------------------------------------------------
# Step 3: Extract Formatted Effect Sizes (SNP:BP:A1:A2, A1, Beta)
# ------------------------------------------------------------------------------
echo "Step 3: Formatting SNP weights for PLINK scoring..."
awk '(FNR>1){print $3":"$4":"$5":"$6,$5,$8"\n"$3":"$4":"$6":"$5,$5,$8}' trait_fullLDM.snpRes > trait_betas.txt

# ------------------------------------------------------------------------------
# Step 4: Individual-level Scoring Across Chromosomes (PLINK 2)
# Note: Designed to run as an HPC job array across chromosomes 1-22
# ------------------------------------------------------------------------------
echo "Step 4: Computing per-chromosome polygenic risk scores..."
for chr in {1..22}; do
    plink_2.00 \
        --vcf ${IMPUTED_VCF_DIR}/chr${chr}.dose.vcf.gz dosage=DS \
        --double-id \
        --extract ${QC_VARIANTS} \
        --score trait_betas.txt cols=+scoresums \
        --threads ${THREADS} \
        --memory ${MEM_MB} \
        --out ${OUTPUT_DIR}/PRS_chr${chr}
done

echo "SBayesRC pipeline completed."
