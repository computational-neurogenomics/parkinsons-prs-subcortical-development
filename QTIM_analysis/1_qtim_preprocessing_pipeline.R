# ==============================================================================
# DATA PREPARATION AND STANDARDISATION
# FOR QTIM COHORT
# ==============================================================================

# 1. LOAD LIBRARIES
library(dplyr)
library(purrr)

# ==============================================================================
# 2. USER CONFIGURATION
# ==============================================================================

# Phenotype data
phenotype_file <- "./data/Phenotypic_Data/QTIM_subcortex.dat"
ct_prs_file    <- "./data/PRS_Data/ClumpThresh_PRS.profile"
sbayes_file    <- "./data/PRS_Data/SBayR_PRS.profile"
pc_file        <- "./data/PCs/GWAS_PrincipalComponentsScores.txt"

# Output target for downstream regression models
output_file    <- "./data/output/qtim_input_for_lme1.csv"

# ==============================================================================
# 3. UPLOADING DATA
# ==============================================================================

# Phenotype data
pheno_data <- read.table(phenotype_file, header = TRUE, sep = " ", stringsAsFactors = FALSE, na.strings = "NA", fill = TRUE)
pheno_data$IID <- as.character(pheno_data$IID)

# PRS for PD in QTIM cohort
ct_df <- read.table(ct_prs_file, header = TRUE)
sbayes_df <- read.table(sbayes_file, header = TRUE)

# PCs data
pc_data <- read.table(pc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "NA")
pc_data <- pc_data %>%
  rename(IID = INDID)

# ==============================================================================
# 4. VARIABLES DECLARATION
# ==============================================================================

phenotypes_list <- c('Accum', 'Amyg', 'Brain_Stem', 'Caud', 
                     'Hippo', 'ICV1', 'Pal', 'Put', 'Thal')#, 'VentralDC')
PCs <- paste0('PC', 1:10)
PRS <- c(paste0('S', 1:8), 'SCORESUM')

# ==============================================================================
# 5. MERGING DATASETS
# ==============================================================================

merged_df <- pheno_data %>%
  full_join(ct_df, by = "IID") %>%
  full_join(sbayes_df, by = "IID") %>%
  full_join(pc_data, by = "IID")

# ==============================================================================
# 6. FORMAT FINAL DATAFRAME
# ==============================================================================

merged_df <- merged_df %>%
  # Retaining 'FID' explicitly for the random intercept in downstream LME models
  select(all_of(c('FID', 'IID', 'SEX', 'AGE', phenotypes_list, PCs, PRS))) %>% 
  na.omit() %>% 
  rename(PDPRS = SCORESUM,
         Accumbens = Accum, Amygdala = Amyg, Brainstem = Brain_Stem, 
         Caudate = Caud, Hippocampus = Hippo, ICV = ICV1,
         Pallidum = Pal, Putamen = Put, Thalamus = Thal, 
         #ventralDC = VentralDC,
         Sex = SEX, Age = AGE)

# ==============================================================================
# 7. STANDARDISATION & OUTLIER REMOVAL
# ==============================================================================

# Overwrite vector with the newly formatted region names
phenotypes_list <- c('Accumbens', 'Amygdala', 'Brainstem', 'Caudate', 
                     'Hippocampus', 'ICV', 'Pallidum', 'Putamen', 'Thalamus')#, 'ventralDC')
PRS <- c(paste0('S', 1:8), 'PDPRS')

# Standardise Subcortical Volumes
for (roi in phenotypes_list) {
  currMean <- mean(merged_df[[roi]], na.rm = TRUE)
  currStd <- sd(merged_df[[roi]], na.rm = TRUE)
  
  # Replace outliers (values beyond 4 standard deviations)
  merged_df[[roi]][merged_df[[roi]] > currMean + 4 * currStd | merged_df[[roi]] < currMean - 4 * currStd] <- NA
  # Standardise
  merged_df[[roi]] <- (merged_df[[roi]] - currMean) / currStd
}

# Standardise Polygenic Risk Scores
for (i in PRS){
  currMean <- mean(merged_df[[i]], na.rm = TRUE)
  currStd <- sd(merged_df[[i]], na.rm = TRUE)
  
  # Replace outliers
  merged_df[[i]][merged_df[[i]] > currMean + 4 * currStd | merged_df[[i]] < currMean - 4 * currStd] <- NA
  # Standardise
  merged_df[[i]] <- (merged_df[[i]] - currMean) / currStd
}

# ==============================================================================
# 8. EXPORT FINAL DATAFRAME
# ==============================================================================

merged_df <- na.omit(merged_df)

# Exporting specifically to link with the regression script
write.csv(merged_df, file = output_file, row.names = FALSE)
print("QTIM data preprocessed.")
