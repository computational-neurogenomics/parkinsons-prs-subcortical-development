# ==============================================================================
# DATA PREPARATION, NEUROCOMBAT HARMONISATION, AND STANDARDISATION
# FOR ABCD COHORT 
# ==============================================================================

# 1. LOAD LIBRARIES
library(dplyr)
library(purrr)
library(neuroCombat)

# ==============================================================================
# 2. USER CONFIGURATION
# ==============================================================================

pheno_data_dir <- "./data/Phenotypic_Data/"
prs_data_file   <- "./data/PRS_Data/PD_ABCD_PRS.profile"


# Output target for downstream regression models
output_file    <- "./data/output/abcd_input_for_ols.csv"
harmonised_output_file    <- "./data/output/abcd_input_for_ols_neurocombat.csv"

# ==============================================================================
# 3. VARIABLES DECLARATION
# ==============================================================================

phenotypes_list <- c('Accumbens', 'Amygdala', 'Brainstem', 'Caudate', 'Hippocampus', 'ICV', 'Pallidum', 'Putamen', 'Thalamus', 'ventralDC')
PCs <- paste0('PC', 1:20)
PRS <- c(paste0('S', 1:8), "SbR")

# ==============================================================================
# 4. READING DATA
# ==============================================================================

prs_df <- read.table(prs_data_file, header = TRUE) # PD PRS Scores in the ABCD cohort

# FORMAT IID COLUMN
prs_df <- prs_df %>%
  mutate(IID = ifelse(grepl("NDAR", IID), sub(".*(NDAR_.*)", "\\1", IID), IID))

data_list <- list() # Initialize an empty list to store the datasets
pheno_files_list <- list.files(path = pheno_data_dir, full.names = TRUE)
for (pheno in pheno_files_list){
  #reading every table in the path, all of them must be tab separated
  pheno_data <- read.table(pheno, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "NA")
  data_list[[pheno]] <- pheno_data #directly accessing the element of data_list that corresponds to the name or index pheno
}

# ==============================================================================
# 5. MERGE DATASETS
# ==============================================================================

pheno_data <- reduce(data_list, function(x, y) full_join(x, y, by = c("IID")))
merged_df <- prs_df %>%
  full_join(pheno_data, by = "IID")

cat("\n---> INITIAL NUMBER OF INDIVIDUALS:", nrow(merged_df), "<---\n\n")

# ==============================================================================
# 6. FORMATTING MERGED DATASET
# ==============================================================================
merged_df <- merged_df %>%
  # we select all variable sof interest + all columns srarting with "ScannerID_"
  select(all_of(c('IID', 'Sex', 'Age', phenotypes_list, PCs, PRS)), starts_with("ScannerID_"))%>% #Select columns of interest
  filter(Sex != 3) %>%
  na.omit() %>% # drop NAs
  rename(PDPRS=SbR)# SBayesR scores

# ==============================================================================
# 7. ComBAT HARMONISATION
# ==============================================================================

# We define which variables to protect
merged_df_neurocombat <-  merged_df
mod <- model.matrix (~Age + Sex + PDPRS, data = merged_df_neurocombat)

# Rebuild the batch vector
scanner_cols <- grep("ScannerID_", names(merged_df_neurocombat), value = TRUE)

# Create a single column looking for which scan has the "1" for each subject. 
# If all of them are 0, the significance of the sibjects belongs to the reference scanner omitted in dummies

merged_df_neurocombat$Single_Scanner_ID <- apply(merged_df_neurocombat[scanner_cols], 1, function(x){
  site <- names(x)[x==1]
  if(length(site) == 0) return("Reference_Scanner")
  return(site[1])
}) 

# extract vector for ComBAT
batch <- merged_df_neurocombat$Single_Scanner_ID

# Create the data matrix
dat <- t(as.matrix(merged_df_neurocombat[, phenotypes_list]))

# Run neuroCombat
print("Running neuroCombat harmonisation...")
harmonised_results <- neuroCombat(dat = dat, batch = batch, mod = mod)

# overwrite raw columns with harmonised data
merged_df_neurocombat[, phenotypes_list] <- t(harmonised_results$dat.combat)
print("Harmonisation complete.")


# ==============================================================================
# 8. STANDARDISATION & OUTLIER REMOVAL
# ==============================================================================

# VARIABLES DECLARATION
phenotypes_list <- c('Accumbens', 'Amygdala', 'Brainstem', 'Caudate', 
                     'Hippocampus', 'ICV', 'Pallidum', 'Putamen', 'Thalamus', 'ventralDC')
PRS <- c(paste0('S', 1:8), 'PDPRS')

# 1. Create a function to handle outlier removal and perfect standardisation
standardise_and_clean <- function(df, columns) {
  for (col in columns) {
    # Calculate initial stats to identify outliers
    currMean <- mean(df[[col]], na.rm = TRUE)
    currStd <- sd(df[[col]], na.rm = TRUE)
    
    # Replace outliers (values beyond ±4 SD) with NA
    df[[col]][df[[col]] > currMean + 4 * currStd | df[[col]] < currMean - 4 * currStd] <- NA
    
    # Recalculate stats
    cleanMean <- mean(df[[col]], na.rm = TRUE)
    cleanStd <- sd(df[[col]], na.rm = TRUE)
    
    # Standardise the column using the clean stats
    df[[col]] <- (df[[col]] - cleanMean) / cleanStd
  }
  return(df)
}

# 2. Apply the function to dataframes
print("Standardising raw data (merged_df)...")
merged_df <- standardise_and_clean(merged_df, c(phenotypes_list, PRS))

print("Standardising harmonised data (merged_df_neurocombat)...")
merged_df_neurocombat <- standardise_and_clean(merged_df_neurocombat, c(phenotypes_list, PRS))

# 3. VERIFY (Checking the harmonised dataframe)
print("Verification of Z-scores for neuroCombat harmonised data:")
sapply(c(phenotypes_list, PRS), function(var) {
  c(mean = mean(merged_df_neurocombat[[var]], na.rm = TRUE),
    sd = sd(merged_df_neurocombat[[var]], na.rm = TRUE))
})

# ==============================================================================
# 9. EXPORT RESULTS
# ==============================================================================

# NON HARMONISED RESULTS 
merged_df <- na.omit(merged_df)
cat("\n---> FINAL NUMBER OF INDIVIDUALS (Non-Harmonised):", nrow(merged_df), "<---\n")
write.csv(merged_df, file = paste0(output_file), row.names = FALSE)

# NeurocomBAT HARMONISED RESULTS 
merged_df_neurocombat <- na.omit(merged_df_neurocombat)
cat("---> FINAL NUMBER OF INDIVIDUALS (Harmonised):", nrow(merged_df_neurocombat), "<---\n\n")
write.csv(merged_df_neurocombat, file = paste0(harmonised_output_file), row.names = FALSE)


