# ==============================================================================
# LINEAR MIXED-EFFECTS MODELS (LME) FOR QTIM COHORT
# Predictors: SBayesRC and Clumping and Thresholding PD PRS
# Outcomes: Subcortical Brain Volumes
# ==============================================================================

# 1. LOAD LIBRARIES
library(lme4)
library(performance)
library(dplyr)
library(broom.mixed)
library(tidyr)
library(purrr)
library(openxlsx)

# ==============================================================================
# 2. USER CONFIGURATION
# ==============================================================================

# Input Data
input_file <- "./data/output/qtim_input_for_lme.csv"

# Output Directory
output_dir <- "./data/output/"

# Core Variables
RoiCodes    <- c("Accumbens","Amygdala","Brainstem","Caudate",
                 "Hippocampus","Pallidum","Putamen","Thalamus",
                 "ventralDC","ICV")
PCcovarStr  <- paste(paste0("PC", 1:10), collapse = "+")
thresholds  <- paste0("S", 1:8)

# ==============================================================================
# 3. UPLOADING DATA
# ==============================================================================

QTIMphenos <- read.csv(input_file)

# ==============================================================================
# 4. FUNCTIONS FOR LINEAR MIXED-EFFECTS MODEL
# ==============================================================================

# A. SBAYESRC FUNCTION
run_sbayes <- function(data, RoiCodes, PCcovarStr, prs_col,
                       model_type = "Model1", adjust_ICV = FALSE) {
  # build base covariate string
  base_cov <- "Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + " 
  
  out <- lapply(RoiCodes, function(roi) {
    # pick reduced vs full
    red_fe  <- paste0(if (adjust_ICV) "ICV + ", base_cov, PCcovarStr)
    full_fe <- paste0(prs_col, " + ", red_fe)
    
    # formulas (FID included as random intercept for family relatedness)
    f_red  <- reformulate(c(strsplit(red_fe, "\\s*\\+\\s*")[[1]], "(1|FID)"), response = roi)
    f_full <- reformulate(c(strsplit(full_fe, "\\s*\\+\\s*")[[1]], "(1|FID)"), response = roi)
    
    # fit
    m_red  <- lmer(f_red,  data = data, REML = FALSE)
    m_full <- lmer(f_full, data = data, REML = FALSE)
    
    # Δ R2 (Marginal)
    r2m <- r2_nakagawa(m_full)$R2_marginal - r2_nakagawa(m_red)$R2_marginal
    
    # extract PRS coefficient
    td <- tidy(m_full, effects="fixed")
    prs_row <- filter(td, term == prs_col)
    if (nrow(prs_row) != 1)
      stop("PRS term not found for ROI ", roi)
    
    tibble(
      ROI       = roi,
      Beta      = prs_row$estimate,
      SE        = prs_row$std.error,
      Tval      = prs_row$statistic,
      P         = 2 * pnorm(-abs(prs_row$statistic)),  # normal approx
      Rsq       = r2m
    )
  }) %>%
    bind_rows()
  
  # print a nice table
  cat("\n", model_type,
      "\nSBayesRC PD PRS -> QTIM subcortical (Δ marginal R²)\n\n")
  out %>%
    mutate(
      `Beta (SE) [t]`            = sprintf("%.2f (%.2f) [%.2f]", Beta, SE, Tval),
      `P`                        = sprintf("%.2e", P),
      `Variance explained (%)`   = sprintf("%.2f", Rsq * 100)
    ) %>%
    select(ROI, `Beta (SE) [t]`, `P`, `Variance explained (%)`) %>%
    print(n = Inf)
  
  invisible(out)
}

# B. CLUMPING & THRESHOLDING FUNCTION 
run_ct <- function(data, RoiCodes, thresholds, PCcovarStr, adjust_ICV = TRUE) {
  base_cov <- "Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + "
  
  results <- list()
  for (thr in thresholds) {
    for (roi in RoiCodes) {
      # build RHS
      red_fe  <- paste0(if (adjust_ICV) "ICV + ", base_cov, PCcovarStr)
      full_fe <- paste0(thr, " + ", if (adjust_ICV) "ICV + ", base_cov, PCcovarStr)
      
      # formulas
      f_red  <- as.formula(paste0(roi, " ~ ", red_fe,  " + (1|FID)"))
      f_full <- as.formula(paste0(roi, " ~ ", full_fe, " + (1|FID)"))
      
      # fit
      m_red  <- lmer(f_red,  data = data, REML = FALSE)
      m_full <- lmer(f_full, data = data, REML = FALSE)
      
      # Δ R² (Marginal)
      R2_red  <- r2_nakagawa(m_red)$R2_marginal
      R2_full <- r2_nakagawa(m_full)$R2_marginal
      deltaR2 <- (R2_full - R2_red) * 100
      
      # extract the threshold term
      td <- tidy(m_full, effects="fixed")
      thr_row <- td[td$term == thr, , drop = FALSE]
      if (nrow(thr_row) != 1) next
      
      results[[length(results)+1]] <- data.frame(
        Threshold             = thr,
        ROI                   = roi,
        Beta                  = thr_row$estimate,
        SE                    = thr_row$std.error,
        Tval                  = thr_row$statistic,
        P                     = 2*pnorm(-abs(thr_row$statistic)),
        Variance_explained_pc = deltaR2,
        stringsAsFactors      = FALSE
      )
    }
  }
  
  df <- bind_rows(results)
  
  # format for printing
  df_print <- df %>%
    mutate(
      `Beta (SE) [t]`          = sprintf("%.2f (%.2f) [%.2f]", Beta, SE, Tval),
      `P-value`                = sprintf("%.2e", P),
      `Variance explained (%)` = sprintf("%.2f", Variance_explained_pc)
    ) %>%
    select(Threshold, ROI, `Beta (SE) [t]`, `P-value`, `Variance explained (%)`)
  
  # simple print
  cat("\nC + T PD PRS predicting QTIM subcortical volumes (Δ marginal R²)\n")
  print(df_print)
  
  invisible(df)
}

# ==============================================================================
# 5. RUNNING MODELS
# ==============================================================================

# RUN SBAYESRC MODELS
results_sbayes_model1 <- run_sbayes(
  data       = QTIMphenos,
  RoiCodes   = RoiCodes,
  PCcovarStr = PCcovarStr,
  prs_col    = "PDPRS",
  model_type = "MODEL 1 (No ICV)",
  adjust_ICV = FALSE
)

results_sbayes_model2 <- run_sbayes(
  data       = QTIMphenos,
  RoiCodes   = RoiCodes,
  PCcovarStr = PCcovarStr,
  prs_col    = "PDPRS",
  model_type = "MODEL 2 (ICV adjusted)",
  adjust_ICV = TRUE
)

# RUN C+T MODELS
ct_model1 <- run_ct(
  data       = QTIMphenos,
  RoiCodes   = RoiCodes,
  thresholds = thresholds,
  PCcovarStr = PCcovarStr,
  adjust_ICV = FALSE
)
ct_model1 <- ct_model1 %>%
  rename(`P-value` = P, `Variance Explained (%)` = Variance_explained_pc)

ct_model2 <- run_ct(
  data       = QTIMphenos,
  RoiCodes   = RoiCodes,
  thresholds = thresholds,
  PCcovarStr = PCcovarStr,
  adjust_ICV = TRUE
)
ct_model2 <- ct_model2 %>%
  rename(`P-value` = P, `Variance Explained (%)` = Variance_explained_pc)

# ==============================================================================
# 6. SENSITIVITY ANALYSIS: Intracranial Volume (ICV)
# ==============================================================================

get_icv_improvement_safe <- function(data, roi, prs_col, PCcovarStr) {
  tryCatch({
    base_cov <- "Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) +"
    
    # 1. El modelo reducido incluye el PRS (sin ICV)
    red_fe  <- paste(prs_col, "+", base_cov, PCcovarStr)
    
    # 2. El modelo completo añade únicamente ICV
    full_fe <- paste(prs_col, "+ ICV +", base_cov, PCcovarStr)
    
    f_red  <- as.formula(paste0(roi, " ~ ", red_fe,  "+ (1|FID)"))
    f_full <- as.formula(paste0(roi, " ~ ", full_fe, "+ (1|FID)"))
    
    m_red  <- lmer(f_red,  data = data, REML = FALSE)
    m_full <- lmer(f_full, data = data, REML = FALSE)
    
    # R2 marginal (Nakagawa)
    r2_red   <- r2_nakagawa(m_red)$R2_marginal
    r2_full  <- r2_nakagawa(m_full)$R2_marginal
    delta_r2 <- r2_full - r2_red
    
    # Prueba de razón de verosimilitud (LRT)
    LRT    <- anova(m_red, m_full)
    chi2   <- LRT$Chisq[2]
    df_val <- if ("Chi Df" %in% names(LRT)) LRT$`Chi Df`[2] else LRT$Df[2]
    pval   <- LRT$`Pr(>Chisq)`[2]
    
    tibble(
      ROI                = roi,
      R2_without_ICV     = r2_red * 100,
      R2_with_ICV        = r2_full * 100,
      Delta_R2           = delta_r2 * 100,
      Chi2               = chi2,
      Degrees_of_freedom = df_val,
      P_value            = pval
    )
  }, error = function(e){
    message(sprintf("Warning: ROI %s failed with error: %s", roi, e$message))
    tibble(
      ROI                = roi,
      R2_without_ICV     = NA_real_,
      R2_with_ICV        = NA_real_,
      Delta_R2           = NA_real_,
      Chi2               = NA_real_,
      Degrees_of_freedom = NA_real_,
      P_value            = NA_real_
    )
  })
}

# Ejecutar únicamente sobre las 9 estructuras subcorticales (excluyendo ICV)
subcortical_rois <- setdiff(RoiCodes, "ICV")
improvements <- map_dfr(subcortical_rois, ~ get_icv_improvement_safe(QTIMphenos, .x, prs_col = "PDPRS", PCcovarStr = PCcovarStr))
print(improvements)

# ==============================================================================
# 7. EXPORT RESULTS
# ==============================================================================

write.csv(results_sbayes_model1, paste0(output_dir, "qtim_results_sbayesrc_model1_rsqmarg.csv"), row.names = FALSE)
write.csv(results_sbayes_model2, paste0(output_dir, "qtim_results_sbayesrc_model2_rsqmarg.csv"), row.names = FALSE)
write.csv(ct_model1, paste0(output_dir, "qtim_results_ct_model1_rsqmarg.csv"), row.names = FALSE)
write.csv(ct_model2, paste0(output_dir, "qtim_results_ct_model2_rsqmarg.csv"), row.names = FALSE)
write.csv(improvements, paste0(output_dir, "qtim_results_icv_sensitivity.csv"), row.names = FALSE)

# Save Master Excel Workbook
wb <- createWorkbook()
addWorksheet(wb, "sbayes_model1")
writeData(wb, "sbayes_model1", results_sbayes_model1)

addWorksheet(wb, "sbayes_model2")
writeData(wb, "sbayes_model2", results_sbayes_model2)

addWorksheet(wb, "CT_model1")
writeData(wb, "CT_model1", ct_model1)

addWorksheet(wb, "CT_model2")
writeData(wb, "CT_model2", ct_model2)

addWorksheet(wb, "ICV_sensitivity")
writeData(wb, "ICV_sensitivity", improvements)

saveWorkbook(wb, file = paste0(output_dir, "QTIM_models_rsqmarg_results.xlsx"), overwrite = TRUE)
print("All QTIM analyses completed.")

