#!/bin/bash
# ==============================================================================
# PIPELINE: CLUMPING & THRESHOLDING (C+T) PRS CALCULATION (PLINK 1.9 / PLINK 2)
# ==============================================================================

REF_PANEL_DIR="/refpanel"
SUMSTATS_INPUT="trait_clumping_input_ids_p.txt"
IMPUTED_PFILE_DIR="/cohort_imputed_plink2_DS"
QC_VARIANTS="/QCed_variants_MAF01_CR90_RSQ60.list"

mkdir -p logfilesClumping PRScalClumping logfilesScoring PRScalScoring

# ------------------------------------------------------------------------------
# Step 1: LD Clumping per Chromosome (r2 < 0.05 within 500 kb)
# ------------------------------------------------------------------------------
echo "Step 1: Running LD Clumping across chromosomes 1-22..."
for chr in {1..22}; do
    plink_2.00 \
        --bim ${REF_PANEL_DIR}/chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.bim_R9names \
        --bed ${REF_PANEL_DIR}/chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.bed \
        --fam ${REF_PANEL_DIR}/chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.ALL.fam \
        --clump ${SUMSTATS_INPUT} \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.05 \
        --clump-kb 500 \
        --threads 1 \
        --memory 60000 \
        --out PRScalClumping/trait_${chr}
done

# ------------------------------------------------------------------------------
# Step 2: Merge Clumped Variant Lists
# ------------------------------------------------------------------------------
echo "Step 2: Merging clumped variants..."
cat PRScalClumping/trait_*.clumps > PRScalClumping/trait_Clumped_All.merged
awk '{if(NR>1) print $3}' PRScalClumping/trait_Clumped_All.merged > clumped_ids.txt
grep -wFf clumped_ids.txt ${SUMSTATS_INPUT} > trait_thresholding_input.txt

# ------------------------------------------------------------------------------
# Step 3: Define Nested P-value Ranges (S1 to S8)
# ------------------------------------------------------------------------------
cat << 'EOF' > p_ranges
S1 0 0.00000005
S2 0 0.00001
S3 0 0.001
S4 0 0.01
S5 0 0.05
S6 0 0.1
S7 0 0.5
S8 0 1.0
EOF

# ------------------------------------------------------------------------------
# Step 4: Multi-threshold Scoring (PLINK 2)
# ------------------------------------------------------------------------------
echo "Step 4: Scoring across p-value thresholds S1-S8..."
for chr in {1..22}; do
    plink_2.00 \
        --pfile ${IMPUTED_PFILE_DIR}_chr${chr} \
        --extract ${QC_VARIANTS} \
        --score trait_thresholding_input.txt 1 2 cols=+scoresums \
        --score-col-nums 3 \
        --q-score-range p_ranges trait_thresholding_input.txt 1 4 \
        --threads 4 \
        --memory 35000 \
        --out PRScalScoring/trait_cohortclumpedPRS_chr${chr}
done

echo "Clumping and thresholding pipeline completed."